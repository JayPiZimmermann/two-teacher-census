"""Publish one replayed auxiliary root after its atlas sweep has exited.

``minimum_locked_aux_root.py`` intentionally refuses publication if the live
sweep's worker accounting changes.  During the final scheduler drain that can
leave a valid, uniquely named ``auxstage`` artifact on disk.  This helper is
the narrow post-sweep counterpart: it requires the named sweep PID and every
local sweep of the requested scope to be gone, takes the ordinary/frontier
root lock, freshly replays the exact stage, and hard-links it into the
canonical name without clobbering anything.

The scheduler restart protocol remains external: do not start another sweep
between this tool's two liveness checks.  A newly visible local sweep makes
the operation fail closed.
"""
import argparse
import fcntl
import glob
import os
import sys

import minimum_atlas_aux_root as A
import minimum_atlas_sweep as S
import minimum_frontier_shard as F
import minimum_locked_aux_root as L


def local_sweep_pids(scope):
    """Return live local atlas sweeps with exactly the requested scope."""
    found = []
    for proc in glob.glob("/proc/[0-9]*"):
        pid = int(os.path.basename(proc))
        stat = A.process_stat(pid)
        argv = A.process_argv(pid)
        index, script = A.local_script(argv, {"minimum_atlas_sweep.py"})
        if stat is None or stat[0] == "Z" or A.process_cwd(pid) != S.HERE \
                or index is None or script != "minimum_atlas_sweep.py":
            continue
        if A.option_value(argv, "--scope", "selected") == scope:
            found.append(pid)
    return sorted(found)


def validate_stage_path(job, path):
    """Require an ordinary auxstage file for exactly ``job`` in this tree."""
    path = os.path.abspath(path)
    if os.path.dirname(path) != os.path.realpath(S.HERE) \
            or os.path.islink(path) or not os.path.isfile(path):
        raise A.PublicationError("stage is not a regular local artifact")
    tag = S.tile_id(job.tile, job.scope)
    suffix = "_%s_r%03d_p%03d.json" % (tag, job.root, job.root)
    if not os.path.basename(path).startswith("minimum_dfs_auxstage_") \
            or not os.path.basename(path).endswith(suffix):
        raise A.PublicationError("stage name does not match the requested root")
    A.validate_staged(job, path)
    return path


def preflight(job, stage, scheduler_pid, sweep_fn=local_sweep_pids,
              pid_exists=os.path.lexists):
    """Validate absence, target isolation, and the stable stage surface."""
    if pid_exists("/proc/%d" % scheduler_pid):
        raise A.AccountingError("named scheduler pid has not exited")
    sweeps = sweep_fn(job.scope)
    if sweeps:
        raise A.AccountingError(
            "local %s sweep still alive: %s" %
            (job.scope, ",".join(map(str, sweeps))))
    tag, canonical, candidate, _, frontier_glob = L.target_paths(job)
    if os.path.lexists(canonical):
        raise A.NoClobberError(
            "canonical shard already exists: %s" % os.path.basename(canonical))
    if os.path.lexists(candidate):
        raise A.NoClobberError(
            "frontier candidate already exists: %s" % os.path.basename(candidate))
    frontier_state = glob.glob(frontier_glob)
    if frontier_state:
        raise A.PublicationError(
            "frontier recovery state already exists: %s" %
            ",".join(sorted(os.path.basename(path)
                            for path in frontier_state)))
    writers = F.live_ordinary_writers(tag, job.root)
    if writers:
        raise A.ActiveTargetError(
            "ordinary shard writer(s) still live for this root: %s" % writers)
    validate_stage_path(job, stage)
    return canonical


def adopt_stage(job, stage, scheduler_pid, checker,
                sweep_fn=local_sweep_pids, replay_fn=L.replay_stage,
                link_fn=os.link, pid_exists=os.path.lexists):
    """Freshly replay and no-clobber publish one post-sweep stage."""
    _, canonical, _, lock_path, _ = L.target_paths(job)
    lock_handle = open(lock_path, "a+")
    try:
        try:
            fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise A.ActiveTargetError(
                "frontier lock is already held: %s" %
                os.path.basename(lock_path)) from exc
        stage = validate_stage_path(job, stage)
        canonical = preflight(job, stage, scheduler_pid, sweep_fn=sweep_fn,
                              pid_exists=pid_exists)
        publication = A.Publication(job, stage, canonical)
        before = L.file_hash(stage)
        digest = replay_fn(publication, checker)
        if digest != before or L.file_hash(stage) != before:
            raise A.PublicationError("stage changed during independent replay")
        # Close normal races from a late writer or cooperative sweep restart.
        preflight(job, stage, scheduler_pid, sweep_fn=sweep_fn,
                  pid_exists=pid_exists)
        stage_stat = os.stat(stage)
        try:
            link_fn(stage, canonical)
        except FileExistsError as exc:
            raise A.NoClobberError(
                "canonical shard appeared during publication: %s" %
                os.path.basename(canonical)) from exc
        try:
            A.validate_staged(job, canonical)
            canonical_stat = os.stat(canonical)
            if (stage_stat.st_dev, stage_stat.st_ino) != \
                    (canonical_stat.st_dev, canonical_stat.st_ino) \
                    or L.file_hash(canonical) != before:
                raise A.PublicationError(
                    "canonical shard differs from replayed stage")
        except BaseException:
            # We alone hold the cooperative root lock.  Roll back only the
            # exact inode linked above; never remove a path replaced by some
            # non-cooperating process after the no-clobber link.
            try:
                current = os.stat(canonical)
                if (current.st_dev, current.st_ino) == \
                        (stage_stat.st_dev, stage_stat.st_ino):
                    os.unlink(canonical)
            except FileNotFoundError:
                pass
            raise
        print("POST-SWEEP PUBLISHED %s <- %s" %
              (os.path.basename(canonical), os.path.basename(stage)),
              flush=True)
        print("POST-SWEEP AUX COMPLETE: sha256 %s" % before, flush=True)
        return canonical, before
    finally:
        lock_handle.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheduler-pid", type=int, required=True)
    parser.add_argument("--scope", choices=("selected", "census"),
                        default="census")
    parser.add_argument("--checker", type=S.parse_dfs_checker,
                        default="minimum_dfs_check_parallel.py")
    parser.add_argument("--stage", required=True)
    parser.add_argument("--root", nargs=5, required=True,
                        metavar=("B0", "B1", "Y0", "Y1", "INDEX"))
    args = parser.parse_args()
    if args.scheduler_pid < 1:
        parser.error("scheduler pid must be positive")
    try:
        job = A.make_job(args.root, args.scope)
        adopt_stage(job, args.stage, args.scheduler_pid, args.checker)
    except Exception as exc:
        print("POST-SWEEP AUX ADOPTION FAILED: %s" % exc,
              file=sys.stderr, flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
