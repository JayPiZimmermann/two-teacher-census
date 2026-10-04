"""Continuously bank pairs of far-queued atlas root shards.

This is deliberately only a scheduler around ``minimum_atlas_aux_root``.  A
batch always contains exactly two distinct tiles and the current first missing
root from each tile.  Each batch is staged under fresh, undiscoverable names
and is published by the guarded STOP/account/no-clobber/CONT transaction in
that module.  The loop stops before the main sweep's frontier gets close, when
the host load reaches the configured ceiling, or after any refused batch.

The sweep journal is neither read nor written.
"""
import argparse
import concurrent.futures
import os
import secrets
import sys
import time

import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S


BATCH_SIZE = 2


def queue_tiles(y_min, beta_min):
    return S.selected_tiles(y_min=y_min, beta_min=beta_min)


def active_frontier(scheduler_pid, expected_workers, scope, queue):
    """Return (active keys, furthest queue index), requiring exact accounting."""
    A.verify_scheduler(scheduler_pid, expected_workers, scope)
    children = A.direct_live_children(scheduler_pid, scope)
    if len(children) != expected_workers:
        raise A.AccountingError(
            "expected %d recognized live children, found %d" %
            (expected_workers, len(children)))
    positions = {S.tile_key(tile): index for index, tile in enumerate(queue)}
    active = {child["tile_key"] for child in children}
    outside = active.difference(positions)
    if outside:
        raise A.AccountingError(
            "live child lies outside the auxiliary queue: %s" %
            (sorted(outside),))
    return active, max(positions[key] for key in active)


def choose_jobs(queue, active, frontier, min_queue_gap, scope):
    """Choose exactly two farthest distinct tiles and their first missing roots."""
    jobs = []
    for index in range(len(queue) - 1, frontier + min_queue_gap, -1):
        tile = queue[index]
        key = S.tile_key(tile)
        if key in active:
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
        jobs.append((index, A.RootJob(key, root, scope)))
        if len(jobs) == BATCH_SIZE:
            return jobs
    return []


def stage_batch(jobs, dfs_script, budget):
    nonce = "%d_%d_%s" % (
        os.getpid(), time.time_ns(), secrets.token_hex(4))
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=BATCH_SIZE) as pool:
        futures = [pool.submit(A.stage_one, job, dfs_script, budget, nonce)
                   for _, job in jobs]
        return [future.result() for future in futures]


def stop_requested(path):
    """A persistent sentinel is sampled only at a completed batch boundary."""
    return path is not None and os.path.lexists(path)


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
    parser.add_argument(
        "--min-queue-gap", type=int, default=12,
        help="stop unless targets are more than this many queue positions "
             "beyond the furthest active main-sweep tile")
    parser.add_argument(
        "--max-load", type=float,
        default=float(os.cpu_count() or 1),
        help="stop before starting a batch at or above this 1-minute load")
    parser.add_argument(
        "--stop-file",
        help="stop cleanly before the next batch when this sentinel exists")
    args = parser.parse_args()
    if args.expected_workers < 1 or args.budget < 1 \
            or args.min_queue_gap < 1 or args.max_load <= 0:
        parser.error("worker count, budget, queue gap, and load must be positive")

    queue = queue_tiles(args.y_min, args.beta_min)
    stop_file = args.stop_file
    if stop_file is not None and not os.path.isabs(stop_file):
        stop_file = os.path.join(S.HERE, stop_file)
    if len(queue) < BATCH_SIZE:
        parser.error("auxiliary queue has fewer than two tiles")
    published = 0
    batches = 0
    try:
        while True:
            if stop_requested(stop_file):
                print("AUXILIARY LOOP STOP: batch-boundary sentinel %s; "
                      "%d shards published in %d batches" %
                      (os.path.basename(stop_file), published, batches),
                      flush=True)
                return 0
            load = os.getloadavg()[0]
            if load >= args.max_load:
                print("AUXILIARY LOOP STOP: load %.2f >= %.2f; "
                      "%d shards published in %d batches" %
                      (load, args.max_load, published, batches), flush=True)
                return 0
            active, frontier = active_frontier(
                args.scheduler_pid, args.expected_workers, args.scope, queue)
            jobs = choose_jobs(queue, active, frontier, args.min_queue_gap,
                               args.scope)
            if len(jobs) != BATCH_SIZE:
                print("AUXILIARY LOOP STOP: fewer than two safe far-queue "
                      "targets; %d shards published in %d batches" %
                      (published, batches), flush=True)
                return 0
            print("AUXILIARY BATCH %d: frontier %d, targets %s, load %.2f" %
                  (batches + 1, frontier,
                   ", ".join("%d:%s:r%d" %
                             (index, S.tile_id(job.tile, args.scope), job.root)
                             for index, job in jobs), load), flush=True)
            publications = stage_batch(jobs, args.dfs_script, args.budget)
            A.guarded_publish(publications, args.scheduler_pid,
                              args.expected_workers, args.scope)
            batches += 1
            published += len(publications)
            print("AUXILIARY BATCH COMPLETE: %d shards total" % published,
                  flush=True)
    except Exception as exc:
        print("AUXILIARY LOOP FAILED after %d shards in %d batches: %s" %
              (published, batches, exc), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
