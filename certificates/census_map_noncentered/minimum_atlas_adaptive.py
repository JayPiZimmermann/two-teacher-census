"""Adaptively certify the F4 selected-minimum atlas.

The fixed atlas has 21 beta columns and 9 y rows.  This driver starts from
coarse unions of those 189 atomic cells.  A coarse rectangle is retained when
its regularity forest and all four physical-boundary forests replay.  If the
regularity DFS exhausts its probe budget, the rectangle is bisected along its
longer grid dimension and the two children are queued.  Thus an unsuccessful
coarse probe is only a scheduling observation, never certificate evidence.

The adaptive journal is separate from ``minimum_atlas_f4_sweep_state.json``.
Already replayed atomic cells in that fixed-grid journal are adopted rather
than recomputed.  The journal and all successful artifacts are restartable.

Typical benchmark run:

  python3 minimum_atlas_adaptive.py --workers 12 --hours 3 --budget 75000

Restartable reflection-source continuation with independently budgeted DFS
shards:

  python3 minimum_atlas_adaptive.py --workers 12 --lower-half --sharded \
      --shard-size 8 --shard-budget 75000 \
      --dfs-script minimum_dfs_cached.py \
      --dfs-checker minimum_dfs_check_cached.py

Omit ``--hours`` to keep refining until the atlas is complete or every
remaining atomic cell has exhausted the probe budget.
"""
import argparse
import concurrent.futures
from decimal import Decimal
import json
import os
import re
import subprocess
import sys
import threading
import time

import minimum_atlas_sweep as fixed


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_PATH = os.path.join(HERE, "minimum_atlas_f4_adaptive_state.json")
LOG_DIR = os.path.join(HERE, "minimum_atlas_f4_adaptive_logs")
FORMAT = "minimum-atlas-adaptive-v1"
PRINT_LOCK = threading.Lock()

BETA = tuple(Decimal("0.1") + Decimal(i) / 10 for i in range(22))
Y = (Decimal("-0.9"), Decimal("-0.85"), Decimal("-0.75"),
     Decimal("-0.65"), Decimal("-0.55"), Decimal("-0.45"),
     Decimal("-0.35"), Decimal("-0.25"), Decimal("-0.15"),
     Decimal("-0.1"))
ROOTS = tuple((bi0, bi1, yi0, yi1)
              for bi0, bi1 in ((0, 7), (7, 14), (14, 21))
              for yi0, yi1 in ((0, 3), (3, 5), (5, 7), (7, 9)))
LOWER_ROOTS = tuple(rect for rect in ROOTS if rect[3] <= 5)
PROGRESS = re.compile(r"nodes\s+(\d+)")


def rect_key(rect):
    return "%d:%d,%d:%d" % rect


def rect_from_key(key):
    fields = re.fullmatch(r"(\d+):(\d+),(\d+):(\d+)", key)
    if fields is None:
        raise ValueError("invalid adaptive rectangle key %s" % key)
    return tuple(map(int, fields.groups()))


def tile_of(rect):
    bi0, bi1, yi0, yi1 = rect
    return {"b0": BETA[bi0], "b1": BETA[bi1],
            "y0": Y[yi0], "y1": Y[yi1]}


def adaptive_id(rect):
    bi0, bi1, yi0, yi1 = rect
    return "f4adapt_bi%02d_%02d_yi%02d_%02d" % (bi0, bi1, yi0, yi1)


def atomic_count(rect):
    bi0, bi1, yi0, yi1 = rect
    return (bi1 - bi0) * (yi1 - yi0)


def split_rect(rect):
    """Bisect on the dimension containing more atomic target cells."""
    bi0, bi1, yi0, yi1 = rect
    nb = bi1 - bi0
    ny = yi1 - yi0
    if nb == ny == 1:
        return ()
    if nb >= ny and nb > 1:
        middle = (bi0 + bi1) // 2
        return ((bi0, middle, yi0, yi1),
                (middle, bi1, yi0, yi1))
    middle = (yi0 + yi1) // 2
    return ((bi0, bi1, yi0, middle),
            (bi0, bi1, middle, yi1))


def load_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def fresh_state():
    return {"format": FORMAT, "strip": fixed.STRIP,
            "passed": {}, "splits": {}, "hard_leaves": {},
            "attempts": []}


def load_state():
    state = load_json(STATE_PATH)
    if state is None:
        return fresh_state()
    if state.get("format") != FORMAT or state.get("strip") != fixed.STRIP:
        raise RuntimeError("adaptive state has an incompatible format/strip")
    return state


def write_state(state):
    temporary = STATE_PATH + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, STATE_PATH)


def fixed_atomic_passes():
    state = load_json(fixed.STATE_PATH) or {}
    passed = state.get("passed", {})
    out = {}
    for bi in range(len(BETA) - 1):
        for yi in range(len(Y) - 1):
            rect = (bi, bi + 1, yi, yi + 1)
            tile = tile_of(rect)
            item = passed.get(fixed.state_key(tile))
            if item is not None:
                out[rect_key(rect)] = item
    return out


def pass_record(rect, item, source):
    tile = tile_of(rect)
    return {"id": item.get("id", adaptive_id(rect)),
            "b0": float(tile["b0"]), "b1": float(tile["b1"]),
            "y0": float(tile["y0"]), "y1": float(tile["y1"]),
            "regularity": item["regularity"],
            "boundaries": item["boundaries"],
            "seconds": item.get("seconds", 0.0), "source": source,
            "atomic_cells": atomic_count(rect)}


def seed_fixed_leaves(state):
    fixed_passed = fixed_atomic_passes()
    for key, item in fixed_passed.items():
        if key not in state["passed"]:
            state["passed"][key] = pass_record(
                rect_from_key(key), item, "fixed-grid-journal")


def frontier(state, roots=ROOTS):
    """Return unresolved leaves of the recorded binary refinement tree."""
    leaves = []

    def visit(rect):
        key = rect_key(rect)
        if key in state["passed"] or key in state["hard_leaves"]:
            return
        children = state["splits"].get(key)
        if children is None:
            leaves.append(rect)
            return
        for child in children:
            visit(tuple(child))

    for root in roots:
        visit(root)
    return leaves


def run_process(argv, timeout, log_path):
    started = time.time()
    try:
        result = subprocess.run(
            argv, cwd=HERE, text=True, stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT, timeout=max(1.0, timeout))
        output = result.stdout
        return result.returncode, output, time.time() - started, False
    except subprocess.TimeoutExpired as exc:
        output = exc.stdout or ""
        if isinstance(output, bytes):
            output = output.decode(errors="replace")
        return 124, output, time.time() - started, True
    finally:
        # A per-attempt log makes progress and failures inspectable without
        # bloating the scheduling journal.
        pass


def save_log(path, argv, code, output, seconds, timed_out):
    os.makedirs(LOG_DIR, exist_ok=True)
    temporary = path + ".tmp"
    with open(temporary, "w") as handle:
        handle.write("command: %s\n" % " ".join(argv))
        handle.write("returncode: %d\nseconds: %.3f\ntimeout: %s\n\n" %
                     (code, seconds, timed_out))
        handle.write(output)
        if output and not output.endswith("\n"):
            handle.write("\n")
    os.replace(temporary, path)


def remaining(deadline):
    return 10 ** 12 if deadline is None else max(0.0, deadline - time.time())


def replay_existing(tile, found_regularity, found_boundaries, deadline,
                    dfs_checker="minimum_dfs_check.py"):
    key = fixed.tile_key(tile)
    regularity = found_regularity.get(key)
    if regularity is None:
        return None
    if remaining(deadline) <= 1:
        raise TimeoutError("benchmark deadline reached before replay")
    fixed.replay_regularity(regularity, dfs_checker)
    boundary_paths = {}
    for face in fixed.FACES:
        path = found_boundaries.get((key, face))
        if path is None:
            return None
        fixed.replay_boundary(path)
        boundary_paths[face] = os.path.basename(path)
    return regularity, boundary_paths


def certify(rect, budget, deadline, found_regularity, found_boundaries,
            sharded=False, shard_size=fixed.SHARD_SIZE,
            shard_budget=fixed.SHARD_NODE_BUDGET, scope="selected",
            dfs_script="minimum_dfs.py",
            dfs_checker="minimum_dfs_check.py"):
    """Attempt one rectangle; only a fully replayed result is successful."""
    started = time.time()
    tile = tile_of(rect)
    key = fixed.tile_key(tile)
    tag = adaptive_id(rect)
    try:
        existing = replay_existing(
            tile, found_regularity, found_boundaries, deadline, dfs_checker)
    except Exception as exc:
        return {"kind": "error", "seconds": time.time() - started,
                "detail": "existing artifact replay failed: %s" % exc}
    if existing is not None:
        regularity, boundary_paths = existing
        return {"kind": "passed", "seconds": time.time() - started,
                "regularity": os.path.basename(regularity),
                "boundaries": boundary_paths, "reused": True}

    if sharded:
        try:
            regularity = fixed.generate_regularity_sharded(
                tile, shard_size=shard_size, shard_budget=shard_budget,
                scope=scope, dfs_script=dfs_script)
            fixed.replay_regularity(regularity, dfs_checker)
        except Exception as exc:
            return {"kind": "error", "seconds": time.time() - started,
                    "detail": "sharded generation/replay failed: %s" % exc}
        observed_nodes = 0
        doc = load_json(regularity)
        if doc is not None:
            observed_nodes = doc.get("stats", {}).get("visited_nodes", 0)
        boundary_paths = {}
        try:
            for face in fixed.FACES:
                if remaining(deadline) <= 1:
                    raise TimeoutError(
                        "benchmark deadline reached on boundaries")
                path = found_boundaries.get((key, face))
                if path is None:
                    path = fixed.generate_boundary(tile, face, scope)
                fixed.replay_boundary(path)
                boundary_paths[face] = os.path.basename(path)
        except Exception as exc:
            return {"kind": "error", "seconds": time.time() - started,
                    "observed_nodes": observed_nodes,
                    "detail": "sharded boundary replay failed: %s" % exc}
        return {"kind": "passed", "seconds": time.time() - started,
                "observed_nodes": observed_nodes,
                "regularity": os.path.basename(regularity),
                "boundaries": boundary_paths, "reused": False,
                "sharded": True}

    argv = fixed.regularity_command(
        tile, tag, 64, 0, 64, budget, scope, dfs_script)
    log_path = os.path.join(LOG_DIR, tag + ".log")
    allowance = remaining(deadline)
    if allowance <= 1:
        return {"kind": "deadline", "seconds": time.time() - started,
                "detail": "benchmark deadline reached before generation"}
    code, output, seconds, timed_out = run_process(argv, allowance, log_path)
    save_log(log_path, argv, code, output, seconds, timed_out)
    progress = [int(match.group(1)) for match in PROGRESS.finditer(output)]
    observed_nodes = max(progress) if progress else 0
    if code != 0:
        kind = "deadline" if timed_out else (
            "budget" if "budget exhausted" in output else "error")
        tail = "\n".join(output.splitlines()[-8:])
        return {"kind": kind, "seconds": time.time() - started,
                "observed_nodes": observed_nodes, "detail": tail}
    try:
        regularity = fixed.artifact_from_output(output)
        fixed.replay_regularity(regularity, dfs_checker)
        boundary_paths = {}
        for face in fixed.FACES:
            if remaining(deadline) <= 1:
                raise TimeoutError("benchmark deadline reached on boundaries")
            path = found_boundaries.get((key, face))
            if path is None:
                path = fixed.generate_boundary(tile, face, scope)
            fixed.replay_boundary(path)
            boundary_paths[face] = os.path.basename(path)
    except Exception as exc:
        return {"kind": "error", "seconds": time.time() - started,
                "observed_nodes": observed_nodes,
                "detail": "post-generation replay/boundary failed: %s" % exc}
    return {"kind": "passed", "seconds": time.time() - started,
            "observed_nodes": observed_nodes,
            "regularity": os.path.basename(regularity),
            "boundaries": boundary_paths, "reused": False}


def record_pass(state, rect, outcome):
    tile = tile_of(rect)
    state["passed"][rect_key(rect)] = {
        "id": adaptive_id(rect),
        "b0": float(tile["b0"]), "b1": float(tile["b1"]),
        "y0": float(tile["y0"]), "y1": float(tile["y1"]),
        "regularity": outcome["regularity"],
        "boundaries": outcome["boundaries"],
        "seconds": outcome["seconds"], "source": "adaptive",
        "atomic_cells": atomic_count(rect)}


def split_or_harden(state, rect, outcome):
    children = split_rect(rect)
    key = rect_key(rect)
    if children:
        state["splits"][key] = [list(child) for child in children]
    else:
        state["hard_leaves"][key] = outcome
    return children


def coverage(state, roots=ROOTS):
    """Count the certified leaves active in the current refinement tree.

    The journal may contain fixed atomic certificates nested strictly inside
    an adaptively successful parent.  They are reusable evidence, but they do
    not constitute additional area and therefore must not be double-counted.
    """
    def visit(rect):
        key = rect_key(rect)
        if key in state["passed"]:
            return atomic_count(rect)
        return sum(visit(tuple(child))
                   for child in state["splits"].get(key, ()))

    return sum(visit(root) for root in roots)


def atlas_boxes(state, roots=ROOTS):
    """The certified leaves of the current refinement tree, as atlas boxes.

    A rectangle that passed is emitted whole and its descendants are not
    visited, exactly as ``coverage`` counts it.  Emitting a nested fixed
    atomic certificate as well would duplicate area and give the overlap
    graph a box strictly inside another.
    """
    boxes = []

    def visit(rect):
        key = rect_key(rect)
        item = state["passed"].get(key)
        if item is not None:
            tile = tile_of(rect)
            boxes.append({
                "id": item.get("id", adaptive_id(rect)),
                "b0": float(tile["b0"]), "b1": float(tile["b1"]),
                "y0": float(tile["y0"]), "y1": float(tile["y1"]),
                "regularity": item["regularity"],
                "boundaries": item["boundaries"]})
            return
        for child in state["splits"].get(key, ()):
            visit(tuple(child))

    for root in roots:
        visit(root)
    return boxes


def target_of(roots):
    return {"b0": float(BETA[min(rect[0] for rect in roots)]),
            "b1": float(BETA[max(rect[1] for rect in roots)]),
            "y0": float(Y[min(rect[2] for rect in roots)]),
            "y1": float(Y[max(rect[3] for rect in roots)])}


def write_atlas(state, path, roots=ROOTS, scope="selected",
                require_complete=True):
    """Assemble the adaptive leaves into an ordinary atlas manifest.

    ``minimum_atlas_sweep.write_atlas`` cannot express this journal: it emits
    one box per fixed atomic tile and demands all of them, while an adaptive
    pass certifies unions of those tiles directly.  The manifest format and
    the authority of ``minimum_atlas_check.py`` are unchanged.
    """
    covered = coverage(state, roots)
    target_cells = sum(atomic_count(rect) for rect in roots)
    if require_complete and covered != target_cells:
        raise RuntimeError("cannot assemble atlas: %d of %d atomic cells are "
                           "certified" % (covered, target_cells))
    boxes = atlas_boxes(state, roots)
    if not boxes:
        raise RuntimeError("cannot assemble atlas: no certified leaf")
    witness = ("minimum_witness_f4probe.json" if scope == "selected"
               else "minimum_witness_f4probe_census.json")
    document = json.load(open(os.path.join(HERE, witness)))
    teacher = document["teacher"]
    base = next((box for box in boxes
                 if box["b0"] <= teacher["beta"] <= box["b1"]
                 and box["y0"] <= teacher["y"] <= box["y1"]), None)
    if base is None:
        raise RuntimeError("no certified box contains the witness teacher")
    manifest = {"format": "minimum-atlas-v1", "scope": scope,
                "strip": fixed.STRIP,
                "witness": witness, "base_box": base["id"], "boxes": boxes}
    if scope == "selected":
        manifest["selected_count"] = document["expected_selected"]
    else:
        manifest["zero_count"] = document["expected_zeros"]
    if require_complete:
        manifest["target"] = target_of(roots)
    destination = path if os.path.isabs(path) else os.path.join(HERE, path)
    temporary = destination + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, destination)
    print("wrote %s with %d certified boxes covering %d/%d atomic cells" %
          (os.path.basename(destination), len(boxes), covered, target_cells))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--hours", type=float)
    parser.add_argument("--budget", type=int, default=75000)
    parser.add_argument("--sharded", action="store_true")
    parser.add_argument("--shard-size", type=int, default=fixed.SHARD_SIZE)
    parser.add_argument("--shard-budget", type=int,
                        default=fixed.SHARD_NODE_BUDGET)
    parser.add_argument(
        "--dfs-script", type=fixed.parse_dfs_script,
        default="minimum_dfs.py",
        help="regularity generator basename; boundary generation and replay "
             "checker selection remain independent")
    parser.add_argument(
        "--dfs-checker", type=fixed.parse_dfs_checker,
        default="minimum_dfs_check.py",
        help="regularity replay checker basename; the cached checker calls "
             "the unchanged ordinary checker after installing the exact "
             "atom cache")
    parser.add_argument(
        "--lower-half", action="store_true",
        help="refine only the 105 source cells with y <= -0.45")
    parser.add_argument("--list", action="store_true")
    parser.add_argument(
        "--scope", default="selected", choices=("selected", "census"),
        help="'selected' counts only zeros that can be local minima; "
             "'census' counts EVERY zero and refuses NEG_T00, which is what "
             "the face-count schema's fold-freeness needs")
    parser.add_argument("--write-atlas")
    parser.add_argument(
        "--partial-atlas", action="store_true",
        help="with --write-atlas, emit the certified leaves so far and omit "
             "the target rectangle, so the manifest claims no coverage")
    args = parser.parse_args()
    if args.workers < 1 or args.budget < 1:
        parser.error("workers and budget must be positive")
    if args.hours is not None and args.hours <= 0:
        parser.error("hours must be positive")
    if args.shard_size < 1 or 64 % args.shard_size \
            or args.shard_budget < 1:
        parser.error("shard-size must divide 64 and shard-budget is positive")

    # Each scope keeps its own journal: the two make different claims about
    # the same rectangles, and a pass under one is not a pass under the other.
    global STATE_PATH, LOG_DIR
    if args.scope != "selected":
        STATE_PATH = os.path.join(
            HERE, "minimum_atlas_f4_adaptive_%s_state.json" % args.scope)
        LOG_DIR = os.path.join(
            HERE, "minimum_atlas_f4_adaptive_%s_logs" % args.scope)
    state = load_state()
    roots = LOWER_ROOTS if args.lower_half else ROOTS
    target_cells = sum(atomic_count(root) for root in roots)
    # --list and --write-atlas are read-only.  Seeding and rewriting the
    # journal here would race a live refinement run, which holds the
    # authoritative state in memory and rewrites the whole file per attempt.
    if args.write_atlas is not None:
        try:
            write_atlas(state, args.write_atlas, roots, args.scope,
                        require_complete=not args.partial_atlas)
        except (RuntimeError, OSError, KeyError) as exc:
            print("ATLAS NOT WRITTEN:", exc)
            return 1
        return 0
    if args.list:
        work = frontier(state, roots)
        print("adaptive frontier: %d rectangles; %d/%d cells covered; "
              "%d hard leaves" %
              (len(work), coverage(state, roots), target_cells,
               len(state["hard_leaves"])))
        for rect in work:
            print(adaptive_id(rect), rect_key(rect), tile_of(rect),
                  "atoms=%d" % atomic_count(rect))
        return 0

    # The fixed-grid journal records SELECTED passes only, so it may not be
    # adopted into a census run.
    if args.scope == "selected":
        seed_fixed_leaves(state)
    write_state(state)
    work = frontier(state, roots)

    deadline = None if args.hours is None else time.time() + args.hours * 3600
    regularity, boundaries = fixed.discover(args.scope)
    queue = list(work)
    running = {}
    completed_attempts = 0
    print("adaptive F4: %d/%d atomic cells already covered, %d coarse "
          "frontier rectangles, %d workers, budget %d, generator %s, "
          "checker %s%s" %
          (coverage(state, roots), target_cells, len(queue), args.workers,
           args.budget, args.dfs_script, args.dfs_checker,
           " for %.2fh" % args.hours if args.hours is not None else ""),
          flush=True)

    with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.workers) as pool:
        while queue or running:
            while queue and len(running) < args.workers \
                    and remaining(deadline) > 1:
                rect = queue.pop(0)
                future = pool.submit(
                    certify, rect, args.budget, deadline,
                    regularity, boundaries, args.sharded,
                    args.shard_size, args.shard_budget, args.scope,
                    args.dfs_script, args.dfs_checker)
                running[future] = rect
                print("START %-35s atoms=%3d beta=[%s,%s] y=[%s,%s]" %
                      (adaptive_id(rect), atomic_count(rect),
                       tile_of(rect)["b0"], tile_of(rect)["b1"],
                       tile_of(rect)["y0"], tile_of(rect)["y1"]),
                      flush=True)
            if not running:
                break
            done, _ = concurrent.futures.wait(
                running, timeout=min(30.0, remaining(deadline) + 2),
                return_when=concurrent.futures.FIRST_COMPLETED)
            if not done:
                if remaining(deadline) <= 0:
                    # Each subprocess has the same deadline timeout and will
                    # return promptly; wait for those outcomes to be journaled.
                    continue
                print("STATUS covered=%d/%d queued=%d running=%d attempts=%d" %
                      (coverage(state, roots), target_cells, len(queue),
                       len(running),
                       completed_attempts), flush=True)
                continue
            for future in done:
                rect = running.pop(future)
                try:
                    outcome = future.result()
                except Exception as exc:
                    outcome = {"kind": "error", "seconds": 0.0,
                               "detail": "worker exception: %s" % exc}
                completed_attempts += 1
                attempt = {"rect": list(rect), "id": adaptive_id(rect),
                           "atomic_cells": atomic_count(rect), **outcome}
                state["attempts"].append(attempt)
                if outcome["kind"] == "passed":
                    record_pass(state, rect, outcome)
                    message = "PASS  %-35s atoms=%3d in %6.0fs coverage=%d/%d" % (
                        adaptive_id(rect), atomic_count(rect),
                        outcome["seconds"], coverage(state, roots),
                        target_cells)
                elif outcome["kind"] == "deadline":
                    # Preserve the unresolved rectangle for a later run.  A
                    # global deadline says nothing about its intrinsic cost.
                    message = "PAUSE %-35s at benchmark deadline after %.0fs" % (
                        adaptive_id(rect), outcome["seconds"])
                else:
                    children = split_or_harden(state, rect, outcome)
                    queue.extend(children)
                    message = "%s %-35s atoms=%3d in %6.0fs -> %s" % (
                        "SPLIT" if children else "HARD ",
                        adaptive_id(rect), atomic_count(rect),
                        outcome["seconds"],
                        ", ".join(adaptive_id(child) for child in children)
                        if children else "atomic leaf")
                write_state(state)
                print(message, flush=True)

    unresolved = frontier(state, roots)
    print("adaptive F4 stopped: %d/%d cells covered, %d unresolved "
          "rectangles, %d hard atomic leaves, %d attempts this run" %
          (coverage(state, roots), target_cells, len(unresolved),
           len(state["hard_leaves"]),
           completed_attempts), flush=True)
    return 0 if coverage(state, roots) == target_cells else 1


if __name__ == "__main__":
    raise SystemExit(main())
