"""Safely precompute disjoint atlas root shards beside a live sweep.

The live sweep snapshots the missing ranges of an active tile only once, so a
second writer must never publish directly into an active tile.  This helper
therefore generates under unique, non-discoverable staging tags.  Publication
is a separate guarded operation:

* stop only the scheduler parent (its numerical children keep running);
* account for exactly the requested number of recognized direct children;
* reject a target tile if any child belongs to it;
* revalidate the exact staged spec, range, scope, and zero UNDECIDED count;
* hard-link into the canonical name with no-clobber semantics; and
* resume the scheduler in a ``finally`` block on every path.

The helper never reads or writes the sweep journal.

Example:
  python3 minimum_atlas_aux_root.py --scheduler-pid 1234 \
      --expected-workers 8 --scope census \
      --root 2.1 2.2 -0.25 -0.15 48 \
      --root 2.0 2.1 -0.15 -0.1 48
"""
import argparse
import concurrent.futures
from dataclasses import dataclass
from decimal import Decimal
import glob
import os
import secrets
import signal
import sys
import threading
import time

import minimum_atlas_sweep as S


PRINT_LOCK = threading.Lock()
CHILD_SCRIPTS = {
    "minimum_dfs.py", "minimum_dfs_cached.py",
    "minimum_dfs_merge.py",
    "minimum_dfs_check.py", "minimum_dfs_check_cached.py",
    "minimum_boundary_dfs.py", "minimum_boundary_check.py",
}


class PublicationError(RuntimeError):
    pass


class ActiveTargetError(PublicationError):
    pass


class NoClobberError(PublicationError):
    pass


class AccountingError(PublicationError):
    pass


@dataclass(frozen=True)
class RootJob:
    tile_key: tuple
    root: int
    scope: str

    @property
    def tile(self):
        return {name: value for name, value in zip(
            ("b0", "b1", "y0", "y1"), self.tile_key)}


@dataclass(frozen=True)
class Publication:
    job: RootJob
    staged_path: str
    canonical_path: str


def process_stat(pid):
    try:
        with open("/proc/%d/stat" % pid) as handle:
            raw = handle.read()
    except OSError:
        return None
    close = raw.rfind(")")
    fields = raw[close + 2:].split() if close >= 0 else []
    if len(fields) < 2:
        return None
    return fields[0], int(fields[1])


def process_argv(pid):
    try:
        with open("/proc/%d/cmdline" % pid, "rb") as handle:
            return [part.decode() for part in handle.read().split(b"\0")
                    if part]
    except (OSError, UnicodeDecodeError):
        return None


def process_cwd(pid):
    try:
        return os.path.realpath(os.readlink("/proc/%d/cwd" % pid))
    except OSError:
        return None


def local_script(argv, allowed):
    if argv is None:
        return None, None
    for index, arg in enumerate(argv):
        name = os.path.basename(arg)
        if name not in allowed:
            continue
        candidate = arg if os.path.isabs(arg) else os.path.join(S.HERE, arg)
        if os.path.realpath(candidate) == os.path.realpath(
                os.path.join(S.HERE, name)):
            return index, name
    return None, None


def option_value(argv, name, default=None):
    try:
        index = argv.index(name)
    except (AttributeError, ValueError):
        return default
    if index + 1 >= len(argv):
        return None
    return argv[index + 1]


def verify_scheduler(pid, expected_workers, scope):
    stat = process_stat(pid)
    argv = process_argv(pid)
    index, name = local_script(argv, {"minimum_atlas_sweep.py"})
    if stat is None or stat[0] == "Z" or process_cwd(pid) != S.HERE \
            or index is None or name != "minimum_atlas_sweep.py":
        raise AccountingError("scheduler pid is not the live local sweep")
    if option_value(argv, "--scope", "selected") != scope:
        raise AccountingError("scheduler scope does not match publication")
    if option_value(argv, "--workers") != str(expected_workers):
        raise AccountingError("scheduler worker count does not match guard")
    if "--sharded" not in argv or option_value(argv, "--shard-size") != "1":
        raise AccountingError("scheduler is not a root-cell sharded sweep")


def wait_stopped(pid, timeout=3.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        stat = process_stat(pid)
        if stat is None or stat[0] == "Z":
            raise AccountingError("scheduler vanished while stopping")
        if stat[0] in ("T", "t"):
            return
        time.sleep(0.02)
    raise AccountingError("scheduler did not enter stopped state")


def tile_lookup():
    return {S.tile_key(tile): tile for tile in S.tiles()}


def checked_tile_key(values, scope):
    try:
        key = S.key_of_bounds(*values)
    except Exception as exc:
        raise AccountingError("child has malformed tile bounds") from exc
    if key not in tile_lookup():
        raise AccountingError("child tile is outside the atlas grid")
    return key


def artifact_tile_key(path, scope):
    if not os.path.isabs(path):
        path = os.path.join(S.HERE, path)
    doc = S.load_json(path)
    spec = doc.get("spec", {}) if isinstance(doc, dict) else {}
    if S.spec_scope(spec) != scope:
        raise AccountingError("checker artifact has the wrong scope")
    return checked_tile_key(
        [spec.get(name) for name in ("b0", "b1", "y0", "y1")], scope)


def recognize_child(pid, scope):
    stat = process_stat(pid)
    argv = process_argv(pid)
    index, script = local_script(argv, CHILD_SCRIPTS)
    if stat is None or stat[0] == "Z" or process_cwd(pid) != S.HERE \
            or index is None:
        raise AccountingError("unrecognized direct child pid %d" % pid)
    args = argv[index + 1:]
    if script in ("minimum_dfs.py", "minimum_dfs_cached.py"):
        if len(args) < 14 or args[12] != "relative" \
                or args[13] != scope:
            raise AccountingError("malformed regularity generator child")
        key = checked_tile_key(args[:4], scope)
    elif script == "minimum_boundary_dfs.py":
        if len(args) < 14 or args[13] != scope:
            raise AccountingError("malformed boundary generator child")
        key = checked_tile_key(args[:4], scope)
    elif script in ("minimum_dfs_check.py",
                    "minimum_dfs_check_cached.py",
                    "minimum_boundary_check.py"):
        if len(args) != 1:
            raise AccountingError("malformed checker child")
        key = artifact_tile_key(args[0], scope)
    elif script == "minimum_dfs_merge.py":
        if len(args) != 1:
            raise AccountingError("malformed merge child")
        matches = [S.tile_key(tile) for tile in S.tiles()
                   if S.tile_id(tile, scope) == args[0]]
        if len(matches) != 1:
            raise AccountingError("merge tag does not name one atlas tile")
        key = matches[0]
    else:  # pragma: no cover - CHILD_SCRIPTS and branches stay in lockstep.
        raise AccountingError("unsupported child script")
    return {"pid": pid, "script": script, "tile_key": key}


def direct_live_children(parent_pid, scope):
    children = []
    for path in glob.glob("/proc/[0-9]*"):
        pid = int(os.path.basename(path))
        stat = process_stat(pid)
        if stat is not None and stat[1] == parent_pid and stat[0] != "Z":
            children.append(recognize_child(pid, scope))
    return sorted(children, key=lambda item: item["pid"])


def validate_staged(job, path):
    doc = S.load_json(path)
    if doc is None or not S.shard_matches(
            path, job.tile, job.root, job.root + 1, job.scope, 64):
        raise PublicationError("staged shard spec/range/scope mismatch")
    stats = doc.get("stats")
    counts = stats.get("counts") if isinstance(stats, dict) else None
    if not isinstance(counts, dict) or counts.get("UNDECIDED") != 0:
        raise PublicationError("staged shard is not complete")
    return doc


def guarded_publish(publications, scheduler_pid, expected_workers, scope):
    if not publications:
        raise PublicationError("nothing to publish")
    targets = [item.job.tile_key for item in publications]
    if len(set(targets)) != len(targets):
        raise PublicationError("auxiliary roots must use disjoint tiles")
    for item in publications:
        validate_staged(item.job, item.staged_path)
    verify_scheduler(scheduler_pid, expected_workers, scope)
    stopped = False
    try:
        os.kill(scheduler_pid, signal.SIGSTOP)
        stopped = True
        wait_stopped(scheduler_pid)
        verify_scheduler(scheduler_pid, expected_workers, scope)
        children = direct_live_children(scheduler_pid, scope)
        if len(children) != expected_workers:
            raise AccountingError(
                "expected %d recognized live children, found %d" %
                (expected_workers, len(children)))
        active = {item["tile_key"] for item in children}
        overlap = active.intersection(targets)
        if overlap:
            raise ActiveTargetError(
                "refusing publication into active tile %s" %
                (sorted(overlap),))
        for item in publications:
            if os.path.lexists(item.canonical_path):
                raise NoClobberError(
                    "canonical shard already exists: %s" %
                    os.path.basename(item.canonical_path))
        for item in publications:
            try:
                os.link(item.staged_path, item.canonical_path)
            except FileExistsError as exc:
                raise NoClobberError(
                    "canonical shard appeared during publication: %s" %
                    os.path.basename(item.canonical_path)) from exc
            validate_staged(item.job, item.canonical_path)
            print("PUBLISHED %s <- %s" %
                  (os.path.basename(item.canonical_path),
                   os.path.basename(item.staged_path)), flush=True)
    finally:
        if stopped and S.process_is_alive(scheduler_pid):
            os.kill(scheduler_pid, signal.SIGCONT)


def make_job(values, scope):
    if len(values) != 5:
        raise ValueError("a root needs b0 b1 y0 y1 index")
    key = tuple(Decimal(value) for value in values[:4])
    if key not in tile_lookup():
        raise ValueError("root tile is outside the atlas grid")
    root = int(values[4])
    if not 0 <= root < 64:
        raise ValueError("root index must lie in [0,64)")
    return RootJob(key, root, scope)


def stage_one(job, dfs_script, budget, nonce):
    tag = "auxstage_%s_%s_r%03d" % (
        nonce, S.tile_id(job.tile, job.scope), job.root)
    expected = S.expected_shard_path(
        job.tile, job.root, job.scope, tag=tag)
    canonical = S.expected_shard_path(job.tile, job.root, job.scope)
    if os.path.lexists(expected):
        raise PublicationError("unique staging path unexpectedly exists")
    with PRINT_LOCK:
        print("STAGING %s root %d as %s" %
              (S.tile_id(job.tile, job.scope), job.root,
               os.path.basename(expected)), flush=True)
    output = S.run_checked(S.regularity_command(
        job.tile, tag, 64, job.root, job.root + 1, budget,
        job.scope, dfs_script))
    emitted = S.artifact_from_output(output)
    if os.path.realpath(emitted) != os.path.realpath(expected):
        raise PublicationError("generator emitted an unexpected staging path")
    doc = validate_staged(job, expected)
    with PRINT_LOCK:
        print("STAGED %s root %d: %d nodes, %.1fs" %
              (S.tile_id(job.tile, job.scope), job.root,
               doc["stats"]["visited_nodes"], doc["stats"]["seconds"]),
              flush=True)
    return Publication(job, expected, canonical)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scheduler-pid", type=int, required=True)
    parser.add_argument("--expected-workers", type=int, default=8)
    parser.add_argument("--scope", choices=("selected", "census"),
                        default="census")
    parser.add_argument("--dfs-script", type=S.parse_dfs_script,
                        default="minimum_dfs_cached.py")
    parser.add_argument("--budget", type=int, default=300000)
    parser.add_argument("--root", action="append", nargs=5, required=True,
                        metavar=("B0", "B1", "Y0", "Y1", "INDEX"))
    args = parser.parse_args()
    if args.expected_workers < 1 or args.budget < 1:
        parser.error("expected-workers and budget must be positive")
    try:
        jobs = [make_job(values, args.scope) for values in args.root]
    except (ValueError, ArithmeticError) as exc:
        parser.error(str(exc))
    if len({job.tile_key for job in jobs}) != len(jobs):
        parser.error("auxiliary roots must use disjoint tiles")
    for job in jobs:
        canonical = S.expected_shard_path(job.tile, job.root, job.scope)
        if os.path.lexists(canonical):
            parser.error("canonical shard already exists: %s" %
                         os.path.basename(canonical))
    nonce = "%d_%d_%s" % (
        os.getpid(), time.time_ns(), secrets.token_hex(4))
    try:
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=len(jobs)) as pool:
            futures = [pool.submit(stage_one, job, args.dfs_script,
                                   args.budget, nonce)
                       for job in jobs]
            publications = [future.result() for future in futures]
        guarded_publish(publications, args.scheduler_pid,
                        args.expected_workers, args.scope)
    except Exception as exc:
        print("AUXILIARY ROOT FAILED: %s" % exc, file=sys.stderr, flush=True)
        return 1
    print("AUXILIARY ROOT COMPLETE: %d shards published" % len(jobs),
          flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
