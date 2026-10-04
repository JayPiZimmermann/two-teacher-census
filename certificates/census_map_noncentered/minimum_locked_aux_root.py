"""Generate one ordinary root shard under the frontier publication lock.

``minimum_atlas_aux_root`` already gives an auxiliary generator a unique,
non-discoverable staging tag and publishes with a guarded, no-clobber hard
link while the live atlas scheduler is stopped and accounted for.  A direct
frontier recovery uses a separate per-root lock.  This wrapper composes both
protocols: it holds that exact frontier lock from before staging until after
guarded publication, and independently replays the staged ordinary shard
before making it canonical.

This is intended for a bounded ordinary-first attempt on one root.  A budget
failure publishes nothing and releases the lock, after which the same root
may be recovered with ``minimum_frontier_shard.py``.
"""
import argparse
import fcntl
import glob
import hashlib
import os
import secrets
import subprocess
import sys
import time

import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S
import minimum_frontier_shard as F


def target_paths(job):
    tag = S.tile_id(job.tile, job.scope)
    canonical = S.expected_shard_path(job.tile, job.root, job.scope)
    lock = os.path.join(
        S.HERE, "minimum_frontier_%s_p%03d.lock" % (tag, job.root))
    candidate = canonical + ".frontier-candidate"
    frontier_glob = os.path.join(
        S.HERE, "minimum_frontier_%s_p%03d_d*" % (tag, job.root))
    return tag, canonical, candidate, lock, frontier_glob


def preflight(job, scheduler_pid, expected_workers, max_load,
              load_one=None):
    """Validate the exact target while its frontier lock is already held."""
    load_one = (lambda: os.getloadavg()[0]) if load_one is None else load_one
    load = load_one()
    if load >= max_load:
        raise A.PublicationError(
            "load %.2f >= %.2f" % (load, max_load))
    tag, canonical, candidate, _, frontier_glob = target_paths(job)
    if os.path.lexists(canonical):
        raise A.NoClobberError(
            "canonical shard already exists: %s" %
            os.path.basename(canonical))
    if os.path.lexists(candidate):
        raise A.NoClobberError(
            "frontier candidate already exists: %s" %
            os.path.basename(candidate))
    frontier_state = glob.glob(frontier_glob)
    if frontier_state:
        raise A.PublicationError(
            "frontier recovery state already exists: %s" %
            ",".join(sorted(os.path.basename(path)
                            for path in frontier_state)))
    writers = F.live_ordinary_writers(tag, job.root)
    if writers:
        raise A.ActiveTargetError(
            "ordinary shard writer(s) already live for this root: %s" %
            writers)
    A.verify_scheduler(scheduler_pid, expected_workers, job.scope)
    children = A.direct_live_children(scheduler_pid, job.scope)
    if len(children) != expected_workers:
        raise A.AccountingError(
            "expected %d recognized live children, found %d" %
            (expected_workers, len(children)))
    if job.tile_key in {child["tile_key"] for child in children}:
        raise A.ActiveTargetError(
            "target tile is active in the live scheduler")
    return load


def file_hash(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def replay_stage(publication, checker):
    """Replay a stable staged shard in a fresh checker process."""
    checker = S.validate_dfs_checker(checker)
    before = file_hash(publication.staged_path)
    result = subprocess.run(
        S.command(checker, publication.staged_path), cwd=S.HERE,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = "\n".join(result.stdout.splitlines()[-30:])
        raise A.PublicationError(
            "independent replay failed (%d):\n%s" %
            (result.returncode, tail))
    after = file_hash(publication.staged_path)
    if after != before:
        raise A.PublicationError("staged shard changed during replay")
    A.validate_staged(publication.job, publication.staged_path)
    print(result.stdout, end="", flush=True)
    print("LOCKED AUX REPLAY GREEN: sha256 %s" % before, flush=True)
    return before


def run_locked(job, scheduler_pid, expected_workers, dfs_script, budget,
               checker, max_load, load_one=None, stage_fn=None,
               replay_fn=None, publish_fn=None):
    """Run the locked stage/replay/publication transaction."""
    stage_fn = A.stage_one if stage_fn is None else stage_fn
    replay_fn = replay_stage if replay_fn is None else replay_fn
    publish_fn = A.guarded_publish if publish_fn is None else publish_fn
    tag, canonical, _, lock_path, _ = target_paths(job)
    lock_handle = open(lock_path, "a+")
    try:
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise A.ActiveTargetError(
                "frontier lock is already held: %s" %
                os.path.basename(lock_path)) from exc
        load = preflight(job, scheduler_pid, expected_workers, max_load,
                         load_one=load_one)
        nonce = "%d_%d_%s" % (
            os.getpid(), time.time_ns(), secrets.token_hex(4))
        print("LOCKED AUX START: %s root %d, load %.2f, budget %d" %
              (tag, job.root, load, budget), flush=True)
        publication = stage_fn(job, dfs_script, budget, nonce)
        if (publication.job != job
                or os.path.realpath(publication.canonical_path)
                != os.path.realpath(canonical)):
            raise A.PublicationError(
                "staging returned the wrong root or canonical path")
        digest = replay_fn(publication, checker)
        publish_fn([publication], scheduler_pid,
                   expected_workers, job.scope)
        if not os.path.lexists(canonical):
            raise A.PublicationError(
                "guarded publication did not create the canonical shard")
        if file_hash(canonical) != digest:
            raise A.PublicationError(
                "canonical shard differs from the replayed stage")
        print("LOCKED AUX COMPLETE: %s sha256 %s" %
              (os.path.basename(canonical), digest), flush=True)
        return publication, digest
    finally:
        lock_handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheduler-pid", type=int, required=True)
    parser.add_argument("--expected-workers", type=int, default=8)
    parser.add_argument("--scope", choices=("selected", "census"),
                        default="census")
    parser.add_argument("--dfs-script", type=S.parse_dfs_script,
                        default="minimum_dfs_cached.py")
    parser.add_argument("--checker", type=S.parse_dfs_checker,
                        default="minimum_dfs_check_parallel.py")
    parser.add_argument("--budget", type=int, default=600000)
    parser.add_argument("--max-load", type=float, default=24.0)
    parser.add_argument("--root", nargs=5, required=True,
                        metavar=("B0", "B1", "Y0", "Y1", "INDEX"))
    args = parser.parse_args()
    if args.expected_workers < 1 or args.budget < 1 or args.max_load <= 0:
        parser.error("worker count, budget, and load ceiling must be positive")
    try:
        job = A.make_job(args.root, args.scope)
        run_locked(job, args.scheduler_pid, args.expected_workers,
                   args.dfs_script, args.budget, args.checker,
                   args.max_load)
    except Exception as exc:
        print("LOCKED AUX FAILED: %s" % exc, file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
