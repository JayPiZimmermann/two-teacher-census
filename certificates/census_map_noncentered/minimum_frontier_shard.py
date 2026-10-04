"""Parallelize one hard minimum-DFS root over its deterministic frontier.

This is a recovery driver for a root cell that is missing because one
``minimum_dfs.py`` invocation exhausted its node budget.  It does not alter
the target spec and therefore reuses all other completed root-cell shards.
Each initial frontier subtree gets an independent process and budget.  If one
still exhausts its budget, only that path is replaced by its two deterministic
children.  This recursive refinement is journaled and resumable.  After every
subtree lands, the driver builds the ordinary
``minimum_dfs_<tag>_pNNN.json`` shard in a temporary file and accepts it only
after a fresh process runs the unchanged checker logic successfully.  The
default replay entry point adds only exact bounded atom memoization; pass
``--checker minimum_dfs_check.py`` for the uncached entry point.
For a large completed tree, ``minimum_dfs_check_parallel.py`` independently
replays a deterministic prefix cover in fresh processes while leaving the
ordinary artifact unchanged.  Its worker count comes from
``MINIMUM_REPLAY_WORKERS`` (default 8).
``minimum_dfs_check_dynamic.py`` preserves those checks but queues smaller
prefix units dynamically so arithmetically expensive units do not pin a
static worker group.

Example (run only after the old writer for p61 has stopped)::

  python3 -u minimum_frontier_shard.py \
    --tag f4_grid_bp050_bp060_ym090_ym085 --cell 61 \
    --depth 8 --workers 16 --budget 600000 --full-width 0.001

The fragment directory is persistent, so interruption and restart lose only
currently running subtrees.
"""
import argparse
import concurrent.futures
import fcntl
import glob
import json
import os
import re
import subprocess
import sys
import time

import minimum_dfs as D
import minimum_dfs_frontier as W


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_FORMAT = "minimum-frontier-recovery-v1"
CHECKERS = ("minimum_dfs_check.py", "minimum_dfs_check_cached.py",
            "minimum_dfs_check_parallel.py",
            "minimum_dfs_check_dynamic.py")


def ordinary_shards(tag):
    return sorted(glob.glob(os.path.join(
        HERE, "minimum_dfs_%s_p*.json" % tag)))


def choose_template(tag, explicit=None):
    candidates = [explicit] if explicit else ordinary_shards(tag)
    for path in candidates:
        if path is None:
            continue
        path = path if os.path.isabs(path) else os.path.join(HERE, path)
        try:
            spec = W.load_spec(path)
        except (OSError, ValueError, json.JSONDecodeError):
            continue
        return path, spec
    raise ValueError("no ordinary shard with an admissible spec for %s" % tag)


def fragment_dir(tag, cell, depth):
    return os.path.join(HERE, "minimum_frontier_%s_p%03d_d%02d" %
                        (tag, cell, depth))


def fragment_path(directory, path):
    return os.path.join(directory, "k%s.json" % (path or "root"))


def state_path(directory):
    return os.path.join(directory, "frontier_state.json")


def load_state(directory, spec, cell, depth):
    path = state_path(directory)
    try:
        with open(path) as handle:
            state = json.load(handle)
    except FileNotFoundError:
        return {"format": STATE_FORMAT, "spec": spec, "cell_idx": cell,
                "base_depth": depth, "split_paths": []}
    if state.get("format") != STATE_FORMAT or state.get("spec") != spec \
            or state.get("cell_idx") != cell \
            or state.get("base_depth") != depth \
            or not isinstance(state.get("split_paths"), list):
        raise ValueError("frontier journal is incompatible with this run")
    for path in state["split_paths"]:
        if not isinstance(path, str) or any(bit not in "01" for bit in path):
            raise ValueError("frontier journal contains an invalid path")
    return state


def write_state(directory, state):
    state = dict(state)
    state["split_paths"] = sorted(set(state["split_paths"]))
    W.write_json_atomic(state_path(directory), state)


def discover_fragments(directory, spec, cell, depth, max_depth):
    """Load structurally complete reusable fragments, ignoring stale files."""
    found = {}
    for path in glob.glob(os.path.join(directory, "k*.json")):
        try:
            with open(path) as handle:
                document = json.load(handle)
            binary_path = document.get("path")
            if not isinstance(binary_path, str) \
                    or not depth <= len(binary_path) <= max_depth \
                    or path != fragment_path(directory, binary_path):
                continue
            fragment = load_fragment(
                path, spec, cell, binary_path)
            if fragment is not None:
                found[binary_path] = fragment
        except (OSError, ValueError, TypeError, json.JSONDecodeError):
            continue
    return found


def plan_frontier(depth, max_depth, available, split_paths):
    """Choose a prefix-free reusable cover and the paths still to run.

    A journaled split takes precedence over a stale parent fragment.  Existing
    descendants cause traversal below an absent parent, so a restart reuses a
    successful sibling and schedules only the missing side.
    """
    if not 0 <= depth <= max_depth:
        raise ValueError("invalid frontier depth range")
    selected = {}
    pending = []
    split_paths = set(split_paths)
    descendant_paths = set(available) | split_paths

    def has_descendant(prefix):
        return any(path.startswith(prefix) and path != prefix
                   for path in descendant_paths)

    def visit(prefix):
        forced_split = prefix in split_paths
        if not forced_split and prefix in available:
            selected[prefix] = available[prefix]
            return
        if len(prefix) < max_depth and (
                forced_split or has_descendant(prefix)):
            visit(prefix + "0")
            visit(prefix + "1")
            return
        pending.append(prefix)

    roots = ([""] if depth == 0 else
             [format(index, "0%db" % depth)
              for index in range(1 << depth)])
    for root in roots:
        visit(root)
    return selected, pending


def load_fragment(path, spec, cell, binary_path):
    try:
        with open(path) as handle:
            fragment = json.load(handle)
        if fragment.get("format") != W.FORMAT \
                or fragment.get("spec") != spec \
                or fragment.get("cell_idx") != cell \
                or fragment.get("path") != binary_path \
                or fragment.get("stats", {}).get("counts", {}).get(
                    "UNDECIDED"):
            return None
        W.audit_fragment(spec, cell, binary_path, fragment)
        return fragment
    except (OSError, ValueError, KeyError, TypeError,
            json.JSONDecodeError):
        return None


def live_ordinary_writers(tag, cell):
    """Find generators that could still publish this ordinary root shard."""
    found = []
    for proc in glob.glob("/proc/[0-9]*/cmdline"):
        try:
            with open(proc, "rb") as handle:
                argv = [part.decode(errors="replace") for part in
                        handle.read().split(b"\0") if part]
            if len(argv) < 2 or os.path.basename(argv[1]) not in (
                    "minimum_dfs.py", "minimum_dfs_cached.py"):
                continue
            index = argv.index(tag)
            lo_idx, hi_idx = int(argv[index + 2]), int(argv[index + 3])
            if lo_idx <= cell < hi_idx:
                found.append(int(proc.split("/")[2]))
        except (OSError, ValueError, IndexError):
            continue
    return sorted(found)


def run_fragment(template, spec, cell, binary_path, budget, full_width,
                 output):
    started = time.time()
    command = [sys.executable, os.path.join(HERE, "minimum_dfs_frontier.py"),
               template, str(cell), binary_path or "-", str(budget), output,
               "--full-width", repr(full_width)]
    result = subprocess.run(command, cwd=HERE, text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    fragment = load_fragment(output, spec, cell, binary_path) \
        if result.returncode == 0 else None
    return binary_path, result.returncode, fragment, result.stdout, \
        time.time() - started


def replay(path, checker="minimum_dfs_check_cached.py"):
    result = subprocess.run(
        [sys.executable, os.path.join(HERE, checker), path],
        cwd=HERE, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError("independent replay failed:\n%s" % result.stdout)
    return result.stdout


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", required=True)
    parser.add_argument("--cell", required=True, type=int)
    parser.add_argument("--depth", type=int, default=8)
    parser.add_argument(
        "--max-depth", type=int, default=20,
        help="recursively split only budget-failed paths up to this depth")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--budget", type=int, default=600000)
    parser.add_argument("--full-width", type=float, default=0.001)
    parser.add_argument("--template")
    parser.add_argument("--checker", choices=CHECKERS,
                        default="minimum_dfs_check_cached.py")
    args = parser.parse_args()
    if args.depth < 0 or args.max_depth < args.depth \
            or args.max_depth > 32 or args.workers < 1 \
            or args.budget < 1 or args.full_width < 0:
        parser.error("require 0 <= depth <= max-depth <= 32, positive "
                     "workers/budget, and nonnegative full-width")
    if re.fullmatch(r"[A-Za-z0-9_]+", args.tag) is None:
        parser.error("tag must contain only letters, digits, and underscores")
    writers = live_ordinary_writers(args.tag, args.cell)
    if writers:
        print("FRONTIER REJECTED: ordinary shard writer(s) still live for "
              "this cell: %s" % writers)
        return 2
    lock_path = os.path.join(
        HERE, "minimum_frontier_%s_p%03d.lock" % (args.tag, args.cell))
    lock_handle = open(lock_path, "a+")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("FRONTIER REJECTED: another frontier driver holds %s" %
              os.path.basename(lock_path))
        return 2
    try:
        template, spec = choose_template(args.tag, args.template)
        if not 0 <= args.cell < spec["n_parts"]:
            raise ValueError("cell is outside the template partition")
    except (OSError, ValueError) as exc:
        print("FRONTIER REJECTED:", exc)
        return 2

    directory = fragment_dir(args.tag, args.cell, args.depth)
    os.makedirs(directory, exist_ok=True)
    try:
        state = load_state(directory, spec, args.cell, args.depth)
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        print("FRONTIER REJECTED:", exc)
        return 2
    available = discover_fragments(
        directory, spec, args.cell, args.depth, args.max_depth)
    fragments, missing_paths = plan_frontier(
        args.depth, args.max_depth, available, state["split_paths"])
    print("frontier p%03d base-depth %d max-depth %d: %d reusable, "
          "%d to run; workers=%d budget=%d full-width=%g" %
          (args.cell, args.depth, args.max_depth, len(fragments),
           len(missing_paths), args.workers, args.budget, args.full_width),
          flush=True)

    failures = []
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.workers) as pool:
        pending = {}

        def submit(binary_path):
            output = fragment_path(directory, binary_path)
            future = pool.submit(
                run_fragment, template, spec, args.cell, binary_path,
                args.budget, args.full_width, output)
            pending[future] = binary_path

        for binary_path in missing_paths:
            submit(binary_path)
        done = 0
        while pending:
            completed, _ = concurrent.futures.wait(
                pending, return_when=concurrent.futures.FIRST_COMPLETED)
            for future in completed:
                pending.pop(future)
                binary_path, code, fragment, output, seconds = future.result()
                done += 1
                if fragment is not None:
                    fragments[binary_path] = fragment
                    print("PASSED path %s nodes=%d %.0fs "
                          "(%d done, %d running)" %
                          (binary_path or "-",
                           fragment["stats"]["visited_nodes"], seconds,
                           done, len(pending)), flush=True)
                    continue
                budget_failed = "budget exhausted" in output
                if budget_failed:
                    state["split_paths"].append(binary_path)
                    write_state(directory, state)
                if budget_failed and len(binary_path) < args.max_depth:
                    submit(binary_path + "0")
                    submit(binary_path + "1")
                    print("SPLIT  path %s after %.0fs -> depth %d "
                          "(%d done, %d running)" %
                          (binary_path or "-", seconds,
                           len(binary_path) + 1, done, len(pending)),
                          flush=True)
                    continue
                failures.append(binary_path)
                tail = "\n".join(output.splitlines()[-3:])
                print("FAILED path %s code=%d %.0fs %s "
                      "(%d done, %d running)" %
                      (binary_path or "-", code, seconds, tail, done,
                       len(pending)), flush=True)
    if failures:
        print("FRONTIER OPEN: %d paths failed: %s" %
              (len(failures), ",".join(failures[:32])))
        return 1

    try:
        merged = W.merge_fragments(
            spec, args.cell, args.depth,
            list(fragments.values()))
        output = os.path.join(HERE, "minimum_dfs_%s_p%03d.json" %
                              (args.tag, args.cell))
        candidate = output + ".frontier-candidate"
        W.write_json_atomic(candidate, merged)
        replay_output = replay(candidate, args.checker)
        os.replace(candidate, output)
    except (OSError, ValueError, RuntimeError) as exc:
        print("FRONTIER MERGE FAILED:", exc)
        return 1
    print(replay_output, end="")
    print("FRONTIER CLOSED: %s (%d bits, %s replayed)" %
          (os.path.basename(output), merged["n_nodes"], args.checker))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
