"""Keep up to two guarded auxiliary root lanes busy without batch skew.

Each completed stage is published serially through
``minimum_atlas_aux_root.guarded_publish``.  Only after that transaction has
resumed the main scheduler may the freed lane be refilled, and every refill
independently requires:

* one-minute host load strictly below the configured ceiling;
* exact recognition of the main scheduler and all of its workers;
* a target beyond the configured queue gap and absent from active workers;
* a tile distinct from the other auxiliary lane; and
* an absent canonical path for that tile's current first missing root.

A persistent stop sentinel is sampled before every refill.  Existing lane
work is allowed to finish and publish, so an intentional later scheduler stop
cannot be undone by this process.  The sweep journal is never read or written.
"""
import argparse
import concurrent.futures
import os
import secrets
import sys
import time

import minimum_atlas_aux_loop as L
import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S


MAX_LANES = 2


def choose_one(queue, excluded, active, frontier, min_queue_gap, scope,
               preferred=None):
    """Return one safe far-queue first-missing root, preferring its old tile."""
    positions = {S.tile_key(tile): index for index, tile in enumerate(queue)}
    order = []
    if preferred is not None and preferred in positions:
        order.append(positions[preferred])
    order.extend(index for index in range(len(queue) - 1, -1, -1)
                 if index not in order)
    for index in order:
        if index <= frontier + min_queue_gap:
            continue
        tile = queue[index]
        key = S.tile_key(tile)
        if key in excluded or key in active:
            continue
        gaps = S.missing_ranges(tile, scope)
        if not gaps:
            continue
        root = gaps[0][0]
        canonical = S.expected_shard_path(tile, root, scope)
        if os.path.lexists(canonical):
            raise A.NoClobberError(
                "first missing root has an occupied canonical path: %s" %
                os.path.basename(canonical))
        return index, A.RootJob(key, root, scope)
    return None


def stage_one(indexed_job, dfs_script, budget):
    index, job = indexed_job
    nonce = "%d_%d_%s" % (
        os.getpid(), time.time_ns(), secrets.token_hex(4))
    return index, A.stage_one(job, dfs_script, budget, nonce)


def run_lanes(queue, scheduler_pid, expected_workers, scope, dfs_script,
              budget, min_queue_gap, max_load, stop_file,
              load_one=None):
    """Run the two-lane coordinator; return (published, started, reason)."""
    load_one = (lambda: os.getloadavg()[0]) if load_one is None else load_one
    published = 0
    started = 0
    reason = "no safe target"
    pending = {}

    def refill(pool, lane, preferred=None):
        nonlocal started, reason
        if L.stop_requested(stop_file):
            reason = "batch-boundary sentinel"
            return False
        load = load_one()
        if load >= max_load:
            reason = "load %.2f >= %.2f" % (load, max_load)
            return False
        active, frontier = L.active_frontier(
            scheduler_pid, expected_workers, scope, queue)
        excluded = {entry[2].tile_key for entry in pending.values()}
        chosen = choose_one(queue, excluded, active, frontier,
                            min_queue_gap, scope, preferred)
        if chosen is None:
            reason = "fewer than one safe far-queue target"
            return False
        index, job = chosen
        print("AUXILIARY LANE %d START: frontier %d, target %d:%s:r%d, "
              "load %.2f" %
              (lane, frontier, index, S.tile_id(job.tile, scope), job.root,
               load), flush=True)
        future = pool.submit(stage_one, chosen, dfs_script, budget)
        pending[future] = (lane, index, job)
        started += 1
        return True

    with concurrent.futures.ThreadPoolExecutor(max_workers=MAX_LANES) as pool:
        for lane in range(MAX_LANES):
            refill(pool, lane)
        while pending:
            done, _ = concurrent.futures.wait(
                tuple(pending), return_when=concurrent.futures.FIRST_COMPLETED)
            for future in done:
                lane, old_index, old_job = pending.pop(future)
                index, publication = future.result()
                if index != old_index or publication.job != old_job:
                    raise A.PublicationError("lane stage returned the wrong job")
                A.guarded_publish(
                    [publication], scheduler_pid, expected_workers, scope)
                published += 1
                print("AUXILIARY LANE %d COMPLETE: %d shards total" %
                      (lane, published), flush=True)
                refill(pool, lane, preferred=old_job.tile_key)
    return published, started, reason


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheduler-pid", type=int, required=True)
    parser.add_argument("--expected-workers", type=int, default=8)
    parser.add_argument("--scope", choices=("selected", "census"),
                        default="census")
    parser.add_argument("--dfs-script", type=S.parse_dfs_script,
                        default="minimum_dfs_cached.py")
    parser.add_argument("--budget", type=int, default=300000)
    parser.add_argument("--beta-min", default="0.5")
    parser.add_argument("--y-min", default="-0.45")
    parser.add_argument("--min-queue-gap", type=int, default=12)
    parser.add_argument("--max-load", type=float, default=24.0)
    parser.add_argument("--stop-file")
    args = parser.parse_args()
    if args.expected_workers < 1 or args.budget < 1 \
            or args.min_queue_gap < 1 or args.max_load <= 0:
        parser.error("worker count, budget, queue gap, and load must be positive")
    stop_file = args.stop_file
    if stop_file is not None and not os.path.isabs(stop_file):
        stop_file = os.path.join(S.HERE, stop_file)
    queue = L.queue_tiles(args.y_min, args.beta_min)
    try:
        published, started, reason = run_lanes(
            queue, args.scheduler_pid, args.expected_workers, args.scope,
            args.dfs_script, args.budget, args.min_queue_gap, args.max_load,
            stop_file)
    except Exception as exc:
        print("AUXILIARY LANES FAILED: %s" % exc, file=sys.stderr, flush=True)
        return 1
    print("AUXILIARY LANES STOP: %s; %d shards published from %d starts" %
          (reason, published, started), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
