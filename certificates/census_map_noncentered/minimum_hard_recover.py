"""Close adaptive hard leaves without retrying known over-budget root cells.

The first escalation pass left two different kinds of missing n=64 roots:

* roots that were never run because the old global queue was stopped;
* roots whose singleton invocation actually exhausted its budget.

Retrying both classes alike wastes workers for many hours.  This driver reads
the old escalation journal but writes a separate recovery journal.  It first
runs only the never-run roots through ``minimum_dfs_cached.py``.  Every
recorded (or newly observed) budget failure is skipped in that phase and then
closed with the recursively refinable ``minimum_frontier_shard.py`` driver.
Completed ordinary shards are the journal of truth, so interruption is safe.

An operator may also intentionally retire an unpublished ordinary root in
favor of exact frontier refinement.  ``--defer-frontier`` records that
scheduling decision separately from an observed budget failure; the two
provenances are never conflated.

Modes::

  python3 minimum_hard_recover.py --list
  python3 minimum_hard_recover.py --list \
      --seed-over-budget f4_grid_bp040_bp050_ym090_ym085#060
  python3 -u minimum_hard_recover.py --normal-only --workers 12
  python3 -u minimum_hard_recover.py --frontier-only --frontier-workers 4
  MINIMUM_REPLAY_WORKERS=8 python3 -u minimum_hard_recover.py \
      --workers 12 --frontier-workers 4 \
      --dfs-checker minimum_dfs_check_parallel.py

The final mode also merges and independently replays any tile whose 64 roots
are complete.  Run only after the old adaptive/escalation schedulers have
exited; the driver refuses while either scheduler is still present.  In the
restart sequence, run ``--normal-only`` BEFORE launching the cached adaptive
sweep.  This ordering is intentional: the startup guard treats even an
adaptive process working on other tiles as a competing future writer.
"""
import argparse
import concurrent.futures
from decimal import Decimal
import glob
import json
import os
import subprocess
import sys
import threading
import time

import minimum_atlas_adaptive as A
import minimum_atlas_sweep as F
import minimum_dfs as D
import minimum_dfs_frontier as W
import minimum_frontier_shard as R
import minimum_shard_escalate as E


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_PATH = os.path.join(HERE, "minimum_hard_recover_state.json")
LOG_DIR = os.path.join(HERE, "minimum_hard_recover_logs")
FORMAT = "minimum-hard-recover-v1"
STATE_LOCK = threading.Lock()
COMPETING_DRIVERS = {"minimum_atlas_adaptive.py",
                     "minimum_shard_escalate.py"}


def fresh_state():
    return {"format": FORMAT, "strip": F.STRIP, "n_parts": 64,
            "cells": {}, "over_budget": {}, "deferred_frontier": {},
            "frontiers": {}, "tiles": {}}


def load_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def load_state():
    state = load_json(STATE_PATH)
    if state is None:
        return fresh_state()
    if state.get("format") != FORMAT or state.get("strip") != F.STRIP \
            or state.get("n_parts") != 64:
        raise ValueError("hard-recovery journal has incompatible metadata")
    # This optional map was added compatibly to the v1 scheduling journal.
    state.setdefault("deferred_frontier", {})
    for key in ("cells", "over_budget", "deferred_frontier", "frontiers",
                "tiles"):
        if not isinstance(state.get(key), dict):
            raise ValueError("hard-recovery journal is missing %s" % key)
    return state


def write_state(state):
    W.write_json_atomic(STATE_PATH, state)


def explicit_tile(values):
    b0, b1, y0, y1 = (Decimal(str(value)) for value in values)
    return {"b0": b0, "b1": b1, "y0": y0, "y1": y1}


def recorded_over_budget(state):
    def normalized(key):
        base, marker, suffix = key.rpartition("#")
        if marker and suffix.isdigit():
            return "%s#%03d" % (base, int(suffix))
        return key

    keys = {normalized(key) for key in state.get("over_budget", {})}
    old = load_json(E.state_path(64)) or {}
    keys.update(normalized(key) for key in old.get("over_budget", {}))
    return keys


def validate_missing_directives(values, tiles, label,
                                missing_fn=E.missing_cells,
                                allow_exact_completed=False,
                                completed_fn=None):
    """Validate scheduling directives without treating them as proof."""
    targets = {E.tile_tag(tile, 64): tile for tile in tiles}
    seeds = set()
    for value in values:
        if not isinstance(value, str) or value.count("#") != 1:
            raise ValueError("%s must have form TAG#CELL" % label)
        tag, suffix = value.split("#")
        if not tag or not suffix.isdigit():
            raise ValueError("%s must have form TAG#CELL" % label)
        cell = int(suffix)
        if not 0 <= cell < 64:
            raise ValueError("%s cell must be in [0,64)" % label)
        tile = targets.get(tag)
        if tile is None:
            raise ValueError("%s is outside the recovery scope: %s" %
                             (label, value))
        missing = set(missing_fn(tile, 64))
        if cell not in missing:
            completed = (exact_completed_root(tile, cell)
                         if completed_fn is None
                         else completed_fn(tile, cell))
            # A direct frontier can publish between a restart wrapper's
            # read-only filter and this validation.  Treat that exact
            # canonical singleton as an idempotent no-op: do not persist its
            # scheduling provenance and do not add it to frontier work.
            if allow_exact_completed and completed:
                continue
            raise ValueError("%s already has a usable shard: %s" %
                             (label, value))
        seeds.add(E.cell_key(tile, cell, 64))
    return seeds


def validate_over_budget_seeds(values, tiles, missing_fn=E.missing_cells,
                               allow_exact_completed=False,
                               completed_fn=None):
    """Validate roots observed to exit on their ordinary node budget."""
    return validate_missing_directives(
        values, tiles, "over-budget seed", missing_fn,
        allow_exact_completed, completed_fn)


def validate_deferred_frontiers(values, tiles, missing_fn=E.missing_cells,
                                allow_exact_completed=False,
                                completed_fn=None):
    """Validate roots intentionally redirected to exact frontier recovery."""
    return validate_missing_directives(
        values, tiles, "deferred frontier", missing_fn,
        allow_exact_completed, completed_fn)


def validate_frontier_directives(seed_values, deferred_values, tiles,
                                 missing_fn=E.missing_cells,
                                 allow_exact_completed=False,
                                 completed_fn=None):
    seeded = validate_over_budget_seeds(
        seed_values, tiles, missing_fn=missing_fn,
        allow_exact_completed=allow_exact_completed,
        completed_fn=completed_fn)
    deferred = validate_deferred_frontiers(
        deferred_values, tiles, missing_fn=missing_fn,
        allow_exact_completed=allow_exact_completed,
        completed_fn=completed_fn)
    overlap = seeded & deferred
    if overlap:
        raise ValueError("roots cannot be both over-budget and deferred: %s" %
                         sorted(overlap))
    return seeded, deferred


def recorded_deferred_frontiers(state):
    """Canonical keys intentionally routed around the ordinary generator."""
    out = set()
    for key in state.get("deferred_frontier", {}):
        base, marker, suffix = key.rpartition("#")
        out.add("%s#%03d" % (base, int(suffix))
                if marker and suffix.isdigit() else key)
    return out


def partition_work(tiles, over_budget, missing_fn=E.missing_cells):
    """Separate missing singleton work from roots already known to fail."""
    normal, frontier = [], []
    for tile in tiles:
        for cell in missing_fn(tile, 64):
            item = (tile, cell)
            if E.cell_key(tile, cell, 64) in over_budget:
                frontier.append(item)
            else:
                normal.append(item)
    return normal, frontier


def competing_driver_pids():
    found = []
    for proc in glob.glob("/proc/[0-9]*/cmdline"):
        try:
            with open(proc, "rb") as handle:
                argv = [part.decode(errors="replace") for part in
                        handle.read().split(b"\0") if part]
            if len(argv) >= 2 and os.path.basename(argv[1]) in \
                    COMPETING_DRIVERS:
                found.append((int(proc.split("/")[2]),
                              os.path.basename(argv[1])))
        except OSError:
            continue
    return sorted(found)


def expected_spec(tile, spec):
    return (D.partition_covers(spec)
            and spec.get("n_parts") == 64
            and spec.get("coord") == "relative"
            and spec.get("scope", "selected") == "selected"
            and all(spec.get(key) == float(tile[key])
                    for key in ("b0", "b1", "y0", "y1"))
            and all(spec.get(key) == F.STRIP[key]
                    for key in ("seam", "delta", "dmax")))


def exact_completed_root(tile, cell):
    """Recognize only the canonical decided singleton for ``tile, cell``."""
    tag = E.tile_tag(tile, 64)
    path = os.path.join(HERE, "minimum_dfs_%s_p%03d.json" % (tag, cell))
    document = load_json(path)
    return (isinstance(document, dict)
            and document.get("lo_idx") == cell
            and document.get("hi_idx") == cell + 1
            and expected_spec(tile, document.get("spec", {}))
            and document.get("stats", {}).get("counts", {}).get(
                "UNDECIDED") == 0)


def save_log(name, command, code, output, seconds):
    os.makedirs(LOG_DIR, exist_ok=True)
    path = os.path.join(LOG_DIR, name)
    document = {"command": command, "returncode": code,
                "seconds": seconds, "output": output}
    W.write_json_atomic(path, document)


def run_normal_cell(tile, cell, budget):
    """Run one never-failed n=64 singleton with the exact-cache wrapper."""
    started = time.time()
    tag = E.tile_tag(tile, 64)
    writers = R.live_ordinary_writers(tag, cell)
    if writers:
        return {"kind": "collision", "seconds": 0.0,
                "detail": "live ordinary writer(s): %s" % writers}
    command = F.command(
        "minimum_dfs_cached.py", tile["b0"], tile["b1"],
        tile["y0"], tile["y1"], F.STRIP["delta"], F.STRIP["dmax"],
        tag, 64, cell, cell + 1, budget, "1e-6", "relative")
    result = subprocess.run(command, cwd=HERE, text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    seconds = time.time() - started
    key = E.cell_key(tile, cell, 64)
    save_log("%s.json" % key.replace("#", "_p"), command,
             result.returncode, result.stdout, seconds)
    if result.returncode:
        return {"kind": ("budget" if "budget exhausted" in result.stdout
                         else "error"),
                "seconds": seconds, "budget": budget,
                "detail": "\n".join(result.stdout.splitlines()[-5:])}
    path = os.path.join(HERE, "minimum_dfs_%s_p%03d.json" % (tag, cell))
    document = load_json(path)
    if document is None or document.get("lo_idx") != cell \
            or document.get("hi_idx") != cell + 1 \
            or not expected_spec(tile, document.get("spec", {})) \
            or document.get("stats", {}).get("counts", {}).get("UNDECIDED"):
        return {"kind": "error", "seconds": seconds,
                "detail": "generator emitted no usable singleton shard"}
    return {"kind": "passed", "seconds": seconds, "budget": budget,
            "nodes": document["stats"].get("visited_nodes", 0),
            "counts": document["stats"].get("counts", {})}


def run_frontier(tile, cell, args):
    tag = E.tile_tag(tile, 64)
    command = F.command(
        "minimum_frontier_shard.py", "--tag", tag, "--cell", cell,
        "--depth", args.depth, "--max-depth", args.max_depth,
        "--workers", args.frontier_workers, "--budget",
        args.frontier_budget,
        "--full-width", args.full_width, "--checker", args.dfs_checker)
    started = time.time()
    result = subprocess.run(command, cwd=HERE, text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    seconds = time.time() - started
    key = E.cell_key(tile, cell, 64)
    save_log("%s_frontier.json" % key.replace("#", "_p"), command,
             result.returncode, result.stdout, seconds)
    remaining = cell in E.missing_cells(tile, 64)
    return {"kind": "passed" if result.returncode == 0 and not remaining
            else "error", "seconds": seconds,
            "detail": "\n".join(result.stdout.splitlines()[-8:])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--frontier-workers", type=int, default=4)
    parser.add_argument("--normal-budget", type=int, default=300000)
    parser.add_argument("--frontier-budget", type=int, default=600000)
    parser.add_argument("--depth", type=int, default=8)
    parser.add_argument("--max-depth", type=int, default=20)
    parser.add_argument("--full-width", type=float, default=0.001)
    parser.add_argument(
        "--dfs-checker", type=F.parse_dfs_checker,
        default="minimum_dfs_check.py",
        help="allowlisted checker used by every frontier publication and "
             "final merged-tile replay; set MINIMUM_REPLAY_WORKERS when "
             "selecting the parallel checker")
    parser.add_argument("--normal-only", action="store_true")
    parser.add_argument("--frontier-only", action="store_true")
    parser.add_argument("--no-finish", action="store_true")
    parser.add_argument("--list", action="store_true")
    parser.add_argument(
        "--seed-over-budget", action="append", default=[],
        metavar="TAG#CELL",
        help="validated missing root observed to exit on budget; repeatable, "
             "and persisted only after the live-scheduler guard passes")
    parser.add_argument(
        "--defer-frontier", action="append", default=[], metavar="TAG#CELL",
        help="validated missing root intentionally retired to exact frontier "
             "recovery; repeatable and distinct from a budget failure")
    parser.add_argument("--tile", nargs=4, action="append",
                        metavar=("B0", "B1", "Y0", "Y1"))
    args = parser.parse_args()
    if args.normal_only and args.frontier_only:
        parser.error("normal-only and frontier-only are mutually exclusive")
    if args.workers < 1 or args.frontier_workers < 1 \
            or args.normal_budget < 1 or args.frontier_budget < 1 \
            or not 0 <= args.depth <= args.max_depth <= 32 \
            or args.full_width < 0:
        parser.error("invalid worker, budget, depth, or width setting")
    try:
        state = load_state()
    except ValueError as exc:
        print("RECOVERY REJECTED:", exc)
        return 2
    tiles = ([explicit_tile(values) for values in args.tile]
             if args.tile else E.hard_leaf_tiles())
    try:
        seeded, deferred_cli = validate_frontier_directives(
            args.seed_over_budget, args.defer_frontier, tiles,
            allow_exact_completed=True)
    except ValueError as exc:
        print("RECOVERY REJECTED:", exc)
        return 2
    over_budget = recorded_over_budget(state) | seeded
    deferred = recorded_deferred_frontiers(state) | deferred_cli
    # An observed budget exit is the more specific provenance if an old
    # hand-edited scheduling journal happened to contain both classifications.
    deferred -= over_budget
    normal, frontier = partition_work(tiles, over_budget | deferred)
    deferred_work = [item for item in frontier
                     if E.cell_key(item[0], item[1], 64) in deferred]
    print("hard recovery: %d tiles, %d missing roots = %d never-failed + "
          "%d frontier (%d over-budget + %d deferred)" %
          (len(tiles), len(normal) + len(frontier), len(normal),
           len(frontier), len(frontier) - len(deferred_work),
           len(deferred_work)), flush=True)
    if args.list:
        for label, work in (("NORMAL", normal), ("FRONTIER", frontier)):
            for tile, cell in work:
                print("%-8s %s" % (label, E.cell_key(tile, cell, 64)))
        return 0
    competitors = competing_driver_pids()
    if competitors:
        print("RECOVERY REJECTED: competing scheduler(s) still live: %s" %
              competitors)
        return 2

    # Persist validated scheduling evidence only after the collision guard.
    # The value is descriptive; ``recorded_over_budget`` deliberately trusts
    # keys only.  Numerical validity still comes solely from later replay.
    if seeded:
        for key in seeded:
            state["over_budget"].setdefault(key, "seeded-cli")
        write_state(state)
    if deferred_cli:
        for key in deferred_cli:
            state["deferred_frontier"].setdefault(key, "intentional-cli")
        write_state(state)

    errors = []
    if not args.frontier_only:
        with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.workers) as pool:
            pending = {pool.submit(run_normal_cell, tile, cell,
                                   args.normal_budget):
                       (tile, cell) for tile, cell in normal}
            done = 0
            for future in concurrent.futures.as_completed(pending):
                tile, cell = pending[future]
                outcome = future.result()
                key = E.cell_key(tile, cell, 64)
                done += 1
                with STATE_LOCK:
                    state["cells"][key] = outcome
                    if outcome["kind"] == "budget":
                        state["over_budget"][key] = args.normal_budget
                    write_state(state)
                print("%-9s %-50s %.0fs (%d/%d)" %
                      (outcome["kind"].upper(), key,
                       outcome.get("seconds", 0), done, len(normal)),
                      flush=True)
                if outcome["kind"] not in ("passed", "budget"):
                    errors.append(key)
    if errors:
        print("RECOVERY OPEN: non-budget errors in %s" % errors[:20])
        return 1
    if args.normal_only:
        return 0

    # Recompute after normal work: successful roots disappear and newly
    # observed budget failures move into the frontier class without retry.
    over_budget = recorded_over_budget(state)
    deferred = recorded_deferred_frontiers(state) - over_budget
    _, frontier = partition_work(tiles, over_budget | deferred)
    for index, (tile, cell) in enumerate(frontier, 1):
        key = E.cell_key(tile, cell, 64)
        outcome = run_frontier(tile, cell, args)
        state["frontiers"][key] = outcome
        write_state(state)
        print("%-9s %-50s %.0fs (%d/%d frontiers)" %
              (outcome["kind"].upper(), key, outcome["seconds"], index,
               len(frontier)), flush=True)
        if outcome["kind"] != "passed":
            errors.append(key)
    if errors:
        print("RECOVERY OPEN: frontier failures in %s" % errors[:20])
        return 1
    if args.no_finish:
        return 0

    for tile in tiles:
        missing = E.missing_cells(tile, 64)
        if missing:
            errors.append("%s missing %s" % (E.tile_tag(tile, 64), missing))
            continue
        outcome = E.finish_tile(tile, 64, args.dfs_checker)
        state["tiles"][F.state_key(tile)] = outcome
        write_state(state)
        print("%-9s %-50s %.0fs" %
              (outcome["kind"].upper(), E.tile_tag(tile, 64),
               outcome["seconds"]), flush=True)
        if outcome["kind"] != "passed":
            errors.append(E.tile_tag(tile, 64))
    print("HARD RECOVERY %s: %d/%d tiles have complete n=64 forests" %
          ("OPEN" if errors else "COMPLETE", len(tiles) - len(errors),
           len(tiles)))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
