"""Finish atomic atlas cells whose grouped root-cell shards exhausted budget.

``minimum_atlas_adaptive.py`` resolves an expensive rectangle by bisecting it
in ``(beta,y)``.  At a 1x1 atomic cell that escape is gone, and the rectangle
is journaled as a hard leaf.  Every hard leaf observed so far fails the same
way: the grouped shard covering root cells ``[32,40)`` or ``[40,48)`` spends
the whole shared budget inside a single root cell.

The shared budget is the only thing that failed.  ``minimum_dfs.py`` takes a
node budget per INVOCATION, not per root cell, so a shard of eight cells and a
shard of one cell are certified identically while the latter gets the whole
budget for its cell.  This driver therefore re-runs only the root-cell ranges
that are still missing from disk, one root cell per invocation, and leaves
every already completed shard in place.

Refinement is exact, not approximate.  ``minimum_dfs_merge.py`` requires the
shard ranges to tile ``[0,n_parts)`` with no gap and no overlap, and shard
files are named by their first root cell, so replacing a grouped shard
``p032 = [32,40)`` by singletons ``p032 = [32,33)``, ``p033``, ..., ``p039``
keeps the partition exact while preserving ``p000`` ... ``p024``.

If one 64-cell root is itself too expensive, ``--n-parts 256`` or ``512``
restarts that TILE under a finer equal root partition.  Its tag, journal, and
log directory include ``n_parts``: shards from different partitions can never
enter one merge.  The resulting full forest is still one ordinary atlas box;
``minimum_dfs_check.py`` derives its partition from the stored spec.

VALIDITY DOMAIN.  This driver schedules and replays; it implements no
numerical verdict and grants no trust.  A cell is recorded as finished only
after ``minimum_dfs_check.py`` replays the merged regularity forest and
``minimum_boundary_check.py`` replays all four physical faces.  A root cell
that exhausts the escalated budget is reported as still open, never as
evidence of anything.  The adaptive and fixed-grid journals are read but
never written, so this is safe to run beside a live sweep.

Usage:
  python3 minimum_shard_escalate.py --list
  python3 minimum_shard_escalate.py --workers 10
  python3 minimum_shard_escalate.py --workers 10 --budget 2000000
  python3 minimum_shard_escalate.py --workers 10 --n-parts 512
  python3 minimum_shard_escalate.py --dfs-script minimum_dfs_cached.py
  python3 minimum_shard_escalate.py --workers 10 --tile 0.1 0.2 -0.9 -0.85

After both this driver and the adaptive driver have stopped, run
``minimum_atlas_reconcile.py --apply`` to adopt fully replayed tiles into the
adaptive journal.  For a refined run pass its separate journal, for example
``minimum_atlas_reconcile.py --escalation-state
minimum_shard_escalate_n512_state.json --apply``.  That command refuses to
race a live adaptive driver.
"""
import argparse
import collections
import concurrent.futures
from decimal import Decimal
import glob
import json
import os
import subprocess
import threading
import time

import minimum_atlas_adaptive as A
import minimum_atlas_sweep as F


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_PATH = os.path.join(HERE, "minimum_shard_escalate_state.json")
LOG_DIR = os.path.join(HERE, "minimum_shard_escalate_logs")
FORMAT = "minimum-shard-escalate-v1"
PRINT_LOCK = threading.Lock()
STATE_LOCK = threading.Lock()
DFS_SCRIPTS = ("minimum_dfs.py", "minimum_dfs_cached.py")


def load_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def state_path(n_parts=64):
    return (STATE_PATH if n_parts == 64 else os.path.join(
        HERE, "minimum_shard_escalate_n%d_state.json" % n_parts))


def log_dir(n_parts=64):
    return (LOG_DIR if n_parts == 64 else os.path.join(
        HERE, "minimum_shard_escalate_n%d_logs" % n_parts))


def fresh_state(n_parts=64):
    return {"format": FORMAT, "strip": F.STRIP, "n_parts": n_parts,
            "cells": {}, "tiles": {}, "over_budget": {}}


def load_state_path(path, expected_n_parts=None):
    state = load_json(path)
    if state is None:
        if expected_n_parts is None:
            raise RuntimeError("escalation state does not exist: %s" % path)
        return fresh_state(expected_n_parts)
    actual_n_parts = state.get("n_parts", 64)
    if state.get("format") != FORMAT or state.get("strip") != F.STRIP:
        raise RuntimeError("escalation state has an incompatible format/strip")
    if expected_n_parts is not None and actual_n_parts != expected_n_parts:
        raise RuntimeError("escalation state has n_parts=%s, expected %s" %
                           (actual_n_parts, expected_n_parts))
    return state


def load_state(n_parts=64):
    return load_state_path(state_path(n_parts), n_parts)


def write_state(state, n_parts=64):
    path = state_path(n_parts)
    temporary = path + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, path)


def tile_tag(tile, n_parts=64):
    base = F.tile_id(tile)
    return base if n_parts == 64 else "%s_n%d" % (base, n_parts)


def missing_cells(tile, n_parts=64):
    """Individual root cells no usable shard covers, in order."""
    return [cell
            for gap_lo, gap_hi in F.missing_ranges(
                tile, n_parts=n_parts, tag=tile_tag(tile, n_parts))
            for cell in range(gap_lo, gap_hi)]


def scheduled_cells(tile, n_parts=64):
    """Missing cells in execution order, without changing their coverage.

    The observed hard region is concentrated at the high-angle end (notably
    old 64-way cell 61).  A refined run has no compatible old shards and must
    eventually compute the whole finer partition, but visiting it backwards
    reaches those inherited hard cells first.  This makes a too-coarse
    refinement fail quickly instead of after hundreds of easy cells.  Keep
    the established 64-way recovery order unchanged.
    """
    cells = missing_cells(tile, n_parts)
    return cells if n_parts == 64 else list(reversed(cells))


def cell_key(tile, cell, n_parts=64):
    return "%s#%03d" % (tile_tag(tile, n_parts), cell)


def generator_script(n_parts=64, override=None):
    """Choose the ordinary generator, except on explicit refined recovery.

    Existing 64-way recovery keeps its established entry point unless the
    caller opts in.  A new finer partition has no compatible work to mix and
    defaults to the bounded exact-cache wrapper.
    """
    if override is not None:
        if override not in DFS_SCRIPTS:
            raise ValueError("unknown DFS script: %s" % override)
        return override
    return "minimum_dfs_cached.py" if n_parts > 64 else "minimum_dfs.py"


def run_cell(tile, cell, budget, n_parts=64, dfs_script=None):
    """Certify exactly one root cell as its own independently budgeted shard."""
    started = time.time()
    tag = tile_tag(tile, n_parts)
    argv = F.command(
        generator_script(n_parts, dfs_script),
        tile["b0"], tile["b1"], tile["y0"], tile["y1"],
        F.STRIP["delta"], F.STRIP["dmax"], tag, n_parts, cell, cell + 1,
        budget, "1e-6", "relative")
    # The node budget, not a wall clock, is what bounds a cell here, so the
    # subprocess runs untimed.  A timeout would report a distinct outcome that
    # says nothing about the cell and would have to be retried anyway.
    result = subprocess.run(argv, cwd=HERE, text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    code, output = result.returncode, result.stdout
    seconds = time.time() - started
    directory = log_dir(n_parts)
    os.makedirs(directory, exist_ok=True)
    A.save_log(os.path.join(directory, "%s_p%03d.log" % (tag, cell)),
               argv, code, output, seconds, False)
    if code != 0:
        kind = "budget" if "budget exhausted" in output else "error"
        return {"kind": kind, "seconds": time.time() - started,
                "budget": budget,
                "detail": "\n".join(output.splitlines()[-4:])}
    path = os.path.join(HERE, "minimum_dfs_%s_p%03d.json" % (tag, cell))
    doc = load_json(path)
    if doc is None or doc.get("lo_idx") != cell \
            or doc.get("hi_idx") != cell + 1:
        return {"kind": "error", "seconds": time.time() - started,
                "budget": budget,
                "detail": "generator did not emit the requested singleton"}
    stats = doc.get("stats", {})
    return {"kind": "passed", "seconds": time.time() - started,
            "budget": budget, "nodes": stats.get("visited_nodes", 0),
            "counts": stats.get("counts", {})}


def merge_tile(tile, n_parts=64):
    """Merge whatever exact shard partition is on disk.

    ``minimum_atlas_sweep.generate_regularity_sharded`` must NOT be used here.
    It expects the uniform eight-cell grouping and would regenerate exactly
    the grouped shard this driver has just replaced by singletons, walking
    back into the shared-budget wall.  The merge script itself imposes no
    grouping: it only demands that the ranges tile ``[0,n_parts)`` exactly.
    """
    tag = tile_tag(tile, n_parts)
    output = F.run_checked(F.command("minimum_dfs_merge.py", tag))
    for line in reversed(output.splitlines()):
        if " -> " in line:
            return os.path.join(HERE, line.split(" -> ", 1)[1].split()[0])
    return os.path.join(HERE, "minimum_dfs_%s_merged.json" % tag)


def finish_tile(tile, n_parts=64, dfs_checker="minimum_dfs_check.py"):
    """Merge, replay, and close the four physical faces of a complete tile."""
    started = time.time()
    try:
        regularity = merge_tile(tile, n_parts)
        F.replay_regularity(regularity, dfs_checker)
    except Exception as exc:
        return {"kind": "error", "seconds": time.time() - started,
                "detail": "merge/replay failed: %s" % exc}
    _, found_boundaries = F.discover()
    key = F.tile_key(tile)
    boundary_paths = {}
    try:
        for face in F.FACES:
            path = found_boundaries.get((key, face))
            if path is None:
                path = F.generate_boundary(tile, face)
            F.replay_boundary(path)
            boundary_paths[face] = os.path.basename(path)
    except Exception as exc:
        return {"kind": "error", "seconds": time.time() - started,
                "regularity": os.path.basename(regularity),
                "detail": "boundary generation/replay failed: %s" % exc}
    return {"kind": "passed", "seconds": time.time() - started,
            "id": tile_tag(tile, n_parts), "n_parts": n_parts,
            "b0": float(tile["b0"]), "b1": float(tile["b1"]),
            "y0": float(tile["y0"]), "y1": float(tile["y1"]),
            "regularity": os.path.basename(regularity),
            "boundaries": boundary_paths}


def hard_leaf_tiles():
    """Atomic tiles the adaptive journal has given up on, read-only."""
    journal = load_json(A.STATE_PATH) or {}
    out = []
    for key in sorted(journal.get("hard_leaves", {})):
        rect = A.rect_from_key(key)
        if A.atomic_count(rect) != 1:
            continue
        out.append(A.tile_of(rect))
    return out


def explicit_tile(values):
    b0, b1, y0, y1 = (Decimal(str(value)) for value in values)
    return {"b0": b0, "b1": b1, "y0": y0, "y1": y1}


def ready_to_finish(tiles, outstanding, n_parts=64):
    """Tiles with a complete on-disk partition and no cell job in flight.

    This includes recovery runs where every singleton landed before an older
    escalation process was interrupted, but the merged tile never reached its
    journal.  Such a tile has no new ``work`` item to trigger the old
    last-cell completion path.
    """
    return [tile for tile in tiles
            if outstanding[tile_tag(tile, n_parts)] == 0
            and not missing_cells(tile, n_parts)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=10)
    parser.add_argument("--budget", type=int, default=F.SHARD_NODE_BUDGET)
    parser.add_argument(
        "--n-parts", type=int, default=64,
        help="equal root partition; a multiple of 64 gives pathological "
             "root cells independent sub-budgets without mixing specs")
    parser.add_argument(
        "--dfs-script", choices=DFS_SCRIPTS,
        help="generator entry point; default minimum_dfs.py at n_parts=64 "
             "and the bounded exact-cache wrapper on refined partitions")
    parser.add_argument("--tile", nargs=4, metavar=("B0", "B1", "Y0", "Y1"),
                        action="append",
                        help="escalate this tile instead of the hard leaves")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    if args.workers < 1 or args.budget < 1 or args.n_parts < 64 \
            or args.n_parts % 64:
        parser.error("workers/budget must be positive and n-parts must be a "
                     "positive multiple of 64")

    tiles = ([explicit_tile(values) for values in args.tile] if args.tile
             else hard_leaf_tiles())
    if not tiles:
        print("no hard atomic leaves to escalate")
        return 0
    # Nearly complete tiles first.  A tile closes only when all root cells
    # are present, so scheduling by remaining work turns the cheapest tiles
    # into closed certificates in minutes instead of parking them behind every
    # expensive cell of every other tile.
    tiles.sort(key=lambda tile: len(missing_cells(tile, args.n_parts)))
    work = [(tile, cell) for tile in tiles
            for cell in scheduled_cells(tile, args.n_parts)]
    if args.list:
        print("escalation targets: %d tiles, %d missing root cells" %
              (len(tiles), len(work)))
        for tile in tiles:
            cells = missing_cells(tile, args.n_parts)
            print("%-42s %2d missing %s" %
                  (tile_tag(tile, args.n_parts), len(cells),
                   "%d..%d" % (cells[0], cells[-1]) if cells else "-"))
        return 0

    state = load_state(args.n_parts)
    dfs_script = generator_script(args.n_parts, args.dfs_script)
    print("escalating %d root cells over %d tiles at n_parts=%d, %d workers, "
          "budget %d, generator %s "
          "(%.1f h per cell at 20 nodes/s)" %
          (len(work), len(tiles), args.n_parts, args.workers, args.budget,
           dfs_script, args.budget / 20.0 / 3600.0), flush=True)
    completed = 0
    closed = 0
    outstanding = collections.Counter(
        tile_tag(tile, args.n_parts) for tile, _ in work)
    finishing = {}
    finisher = concurrent.futures.ThreadPoolExecutor(max_workers=2)
    for tile in ready_to_finish(tiles, outstanding, args.n_parts):
        finishing[finisher.submit(finish_tile, tile, args.n_parts)] = tile
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.workers) as pool:
        pending = {pool.submit(run_cell, tile, cell, args.budget,
                               args.n_parts, dfs_script):
                   (tile, cell) for tile, cell in work}
        for future in concurrent.futures.as_completed(pending):
            tile, cell = pending[future]
            try:
                outcome = future.result()
            except Exception as exc:
                outcome = {"kind": "error", "seconds": 0.0,
                           "detail": "worker exception: %s" % exc}
            completed += 1
            key = cell_key(tile, cell, args.n_parts)
            with STATE_LOCK:
                state["cells"][key] = outcome
                if outcome["kind"] == "budget":
                    state["over_budget"][key] = outcome["budget"]
                else:
                    state["over_budget"].pop(key, None)
                write_state(state, args.n_parts)
            note = ("nodes=%d" % outcome["nodes"]
                    if outcome["kind"] == "passed"
                    else outcome.get("detail", "").strip().splitlines()[-1:])
            with PRINT_LOCK:
                print("%-6s %-46s %6.0fs %s (%d/%d)" %
                      (outcome["kind"].upper(), key, outcome["seconds"],
                       note if isinstance(note, str) else
                       (note[0] if note else ""),
                       completed, len(work)), flush=True)
            # Close a tile the moment its last root cell lands, rather than
            # after the whole queue drains: a tile needing eight cheap cells
            # would otherwise wait behind every expensive cell of every other
            # tile before becoming a certificate.
            tag = tile_tag(tile, args.n_parts)
            outstanding[tag] -= 1
            if outstanding[tag] == 0 \
                    and not missing_cells(tile, args.n_parts):
                finishing[finisher.submit(
                    finish_tile, tile, args.n_parts)] = tile

    for future in concurrent.futures.as_completed(finishing):
        tile = finishing[future]
        try:
            outcome = future.result()
        except Exception as exc:
            outcome = {"kind": "error", "seconds": 0.0,
                       "detail": "finisher exception: %s" % exc}
        with STATE_LOCK:
            state_key = F.state_key(tile)
            if args.n_parts != 64:
                state_key += "#n%d" % args.n_parts
            state["tiles"][state_key] = outcome
            write_state(state, args.n_parts)
        if outcome["kind"] == "passed":
            closed += 1
            print("CLOSED %-42s in %.0fs, regularity and four faces replayed"
                  % (tile_tag(tile, args.n_parts), outcome["seconds"]),
                  flush=True)
        else:
            print("FAILED %-42s %s" %
                  (tile_tag(tile, args.n_parts), outcome.get("detail")),
                  flush=True)
    finisher.shutdown()
    for tile in tiles:
        remaining = missing_cells(tile, args.n_parts)
        if remaining:
            print("OPEN   %-42s %d root cells still over budget: %s" %
                  (tile_tag(tile, args.n_parts), len(remaining), remaining),
                  flush=True)
    print("escalation finished: %d/%d tiles closed" % (closed, len(tiles)),
          flush=True)
    return 0 if closed == len(tiles) else 1


if __name__ == "__main__":
    raise SystemExit(main())
