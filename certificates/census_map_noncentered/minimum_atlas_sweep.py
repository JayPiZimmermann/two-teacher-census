"""Run and replay the rectangular F4 selected-minimum atlas to completion.

The individual certificate programs intentionally emit one small, independently
replayable artifact.  This driver supplies only scheduling and restartability:

* the target rectangle is tiled by exact decimal endpoints;
* existing artifacts are identified from their stored specifications, not from
  filenames;
* a tile is recorded as complete only after its regularity forest and all four
  physical-boundary forests replay successfully;
* failures do not cancel unrelated tiles, and a later run retries them;
* the state file is only a work journal.  ``minimum_atlas_check.py`` remains
  the final authority for an assembled atlas.

No numerical verdict is implemented here.

Usage:
  python3 minimum_atlas_sweep.py --list
  python3 minimum_atlas_sweep.py --workers 20
  python3 minimum_atlas_sweep.py --workers 12 --sharded --lower-half
  python3 minimum_atlas_sweep.py --sharded \
      --dfs-script minimum_dfs_cached.py \
      --dfs-checker minimum_dfs_check_cached.py
  python3 minimum_atlas_sweep.py --workers 12 --adopt-running
  python3 minimum_atlas_sweep.py --workers 20 --limit 20
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


HERE = os.path.dirname(os.path.abspath(__file__))
STATE_PATH = os.path.join(HERE, "minimum_atlas_f4_sweep_state.json")
FACES = ("seam_lo", "seam_hi", "gap_lo", "gap_hi")
STRIP = {"seam": 0.137, "delta": 0.001, "dmax": 3.13}
REGULARITY_NODE_BUDGET = 2000000
SHARD_NODE_BUDGET = 500000
SHARD_SIZE = 8
PRINT_LOCK = threading.Lock()
DFS_SCRIPTS = ("minimum_dfs.py", "minimum_dfs_cached.py")
DFS_CHECKERS = ("minimum_dfs_check.py", "minimum_dfs_check_cached.py",
                "minimum_dfs_check_parallel.py",
                "minimum_dfs_check_dynamic.py")


def decimal_range(first, stop, step):
    value = Decimal(first)
    end = Decimal(stop)
    stride = Decimal(step)
    while value < end:
        yield value
        value += stride


def key_of_bounds(b0, b1, y0, y1):
    return tuple(Decimal(str(value)) for value in (b0, b1, y0, y1))


def tile_key(tile):
    return tile["b0"], tile["b1"], tile["y0"], tile["y1"]


def tile_id(tile, scope="selected"):
    """Artifact tag for a tile.

    The scope is part of the tag because the two claims are different and
    must not share a filename.  Note the asymmetry: a census forest admits
    no NEG_T00 leaf, so it is also a valid SELECTED certificate, while a
    selected forest is not a valid census one.  Keeping them in separate
    files is what stops the weaker artifact being consumed as the stronger.
    """
    if tile_key(tile) == (Decimal("0.7"), Decimal("0.8"),
                          Decimal("-0.45"), Decimal("-0.35")):
        base = "f4probe"
    else:
        def code(value):
            sign = "m" if value < 0 else "p"
            scaled = abs(value * 100)
            return "%s%03d" % (sign, int(scaled))

        base = "f4_grid_b%s_b%s_y%s_y%s" % (
            code(tile["b0"]), code(tile["b1"]),
            code(tile["y0"]), code(tile["y1"]))
    return base if scope == "selected" else "%s_%s" % (base, scope)


def spec_scope(spec):
    """The scope a stored spec claims; absent means the original "selected"."""
    return spec.get("scope", "selected")


def tiles():
    # The vertical grid is anchored at the committed witness box
    # ``[-0.45,-0.35]``.  Two half-height edge rows close the continuation
    # rectangle exactly without throwing away the already replayed frontier.
    y_intervals = [(Decimal("-0.9"), Decimal("-0.85"))]
    y_intervals.extend(
        (y0, y0 + Decimal("0.1"))
        for y0 in decimal_range("-0.85", "-0.15", "0.1"))
    y_intervals.append((Decimal("-0.15"), Decimal("-0.1")))
    out = []
    for b0 in decimal_range("0.1", "2.2", "0.1"):
        b1 = b0 + Decimal("0.1")
        for y0, y1 in y_intervals:
            out.append({"b0": b0, "b1": b1, "y0": y0, "y1": y1})

    def distance_from_base(tile):
        b0, b1, y0, y1 = tile_key(tile)
        if b1 < Decimal("0.7"):
            db = Decimal("0.7") - b1
        elif Decimal("0.8") < b0:
            db = b0 - Decimal("0.8")
        else:
            db = Decimal(0)
        if y1 < Decimal("-0.45"):
            dy = Decimal("-0.45") - y1
        elif Decimal("-0.35") < y0:
            dy = y0 - Decimal("-0.35")
        else:
            dy = Decimal(0)
        return db + dy, b0, y0

    # Grow in connected Manhattan layers.  The ordering is useful on an
    # interrupted run even though workers within a layer finish asynchronously.
    return sorted(out, key=distance_from_base)


def load_json(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def same_strip(spec):
    return all(spec.get(name) == value for name, value in STRIP.items())


def full_forest(doc):
    spec = doc.get("spec", {})
    return (doc.get("lo_idx") == 0
            and doc.get("hi_idx") == spec.get("n_parts"))


def discover(scope="selected"):
    """Full-period forests on disk, restricted to ONE scope.

    The filter matters in one direction only.  A census forest admits no
    NEG_T00 leaf, so it also replays as a selected certificate; a selected
    forest does NOT replay as a census one.  Without the scope test a census
    consumer would happily accept the weaker artifact.
    """
    regularity = {}
    boundaries = {}
    # Full forests may be emitted directly as ``..._p000.json`` or assembled
    # from exact root-cell shards as ``..._merged.json``.  ``full_forest``
    # below rejects the individual shard files.
    for path in glob.glob(os.path.join(HERE, "minimum_dfs_f4*.json")):
        doc = load_json(path)
        if doc is None:
            continue
        spec = doc.get("spec", {})
        if (spec.get("coord") != "relative" or not same_strip(spec)
                or spec_scope(spec) != scope or not full_forest(doc)):
            continue
        key = key_of_bounds(*(spec.get(name) for name in
                              ("b0", "b1", "y0", "y1")))
        regularity[key] = path
    for path in glob.glob(
            os.path.join(HERE, "minimum_boundary_f4*_p000.json")):
        doc = load_json(path)
        if doc is None:
            continue
        spec = doc.get("spec", {})
        face = spec.get("face")
        if (spec.get("coord") != "absolute" or face not in FACES
                or not same_strip(spec) or spec_scope(spec) != scope
                or not full_forest(doc)):
            continue
        key = key_of_bounds(*(spec.get(name) for name in
                              ("b0", "b1", "y0", "y1")))
        boundaries[(key, face)] = path
    return regularity, boundaries


def process_is_alive(pid):
    """Return false for vanished and zombie processes."""
    try:
        with open("/proc/%d/stat" % pid) as handle:
            stat = handle.read()
    except OSError:
        return False
    close = stat.rfind(")")
    return close >= 0 and stat[close + 2:close + 3] != "Z"


def tile_from_generator_process(pid, by_key):
    """Recognize one live allowed relative DFS generator process in HERE."""
    try:
        cwd = os.path.realpath(os.readlink("/proc/%d/cwd" % pid))
        with open("/proc/%d/cmdline" % pid, "rb") as handle:
            argv = [part.decode() for part in handle.read().split(b"\0")
                    if part]
    except (OSError, UnicodeDecodeError):
        return None
    if cwd != os.path.realpath(HERE) or not process_is_alive(pid):
        return None
    def allowed_script(index):
        arg = argv[index]
        name = os.path.basename(arg)
        if name not in DFS_SCRIPTS:
            return False
        candidate = arg if os.path.isabs(arg) else os.path.join(cwd, arg)
        return os.path.realpath(candidate) == os.path.realpath(
            os.path.join(HERE, name))

    script_index = next(
        (index for index in range(len(argv)) if allowed_script(index)), None)
    # coord is the thirteenth argument after the script.  Scope is optional
    # in old selected invocations and follows coord in current invocations,
    # so testing argv[-1] would reject every explicit-scope process.
    if script_index is None or len(argv) <= script_index + 13 \
            or argv[script_index + 13] != "relative":
        return None
    try:
        key = key_of_bounds(*argv[script_index + 1:script_index + 5])
    except Exception:
        return None
    return by_key.get(key)


def discover_running_generators(all_tiles):
    """Map live in-directory F4 generators to their exact atlas tiles."""
    by_key = {tile_key(tile): tile for tile in all_tiles}
    adopted = {}
    for proc in glob.glob("/proc/[0-9]*"):
        pid = int(os.path.basename(proc))
        tile = tile_from_generator_process(pid, by_key)
        if tile is None:
            continue
        key = tile_key(tile)
        if key in adopted:
            raise RuntimeError("multiple live generators for tile %s" %
                               state_key(tile))
        adopted[key] = (pid, tile)
    return adopted


def command(script, *args):
    return [sys.executable, os.path.join(HERE, script)] + [str(x) for x in args]


def validate_dfs_script(script):
    """Accept only the two in-directory regularity generator basenames."""
    if not isinstance(script, str) or os.path.basename(script) != script \
            or script not in DFS_SCRIPTS:
        raise ValueError("DFS script must be one of: %s" %
                         ", ".join(DFS_SCRIPTS))
    return script


def parse_dfs_script(script):
    try:
        return validate_dfs_script(script)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc))


def validate_dfs_checker(script):
    """Accept only explicitly allowlisted in-directory checker basenames."""
    if not isinstance(script, str) or os.path.basename(script) != script \
            or script not in DFS_CHECKERS:
        raise ValueError("DFS checker must be one of: %s" %
                         ", ".join(DFS_CHECKERS))
    return script


def parse_dfs_checker(script):
    try:
        return validate_dfs_checker(script)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc))


def regularity_command(tile, tag, n_parts, lo_idx, hi_idx, budget,
                       scope="selected", dfs_script="minimum_dfs.py"):
    """Construct one regularity command from a validated local entry point."""
    return command(
        validate_dfs_script(dfs_script),
        tile["b0"], tile["b1"], tile["y0"], tile["y1"],
        STRIP["delta"], STRIP["dmax"], tag, n_parts, lo_idx, hi_idx,
        budget, "1e-6", "relative", scope)


def run_checked(argv):
    result = subprocess.run(argv, cwd=HERE, text=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode:
        tail = "\n".join(result.stdout.splitlines()[-20:])
        raise RuntimeError("command failed (%d): %s\n%s" %
                           (result.returncode, " ".join(argv), tail))
    return result.stdout


def artifact_from_output(output):
    for line in reversed(output.splitlines()):
        if line.startswith("wrote "):
            name = line.split()[1]
            return os.path.join(HERE, name)
    raise RuntimeError("generator did not report an artifact")


def generate_regularity(tile, scope="selected",
                        dfs_script="minimum_dfs.py"):
    tag = tile_id(tile, scope)
    output = run_checked(regularity_command(
        tile, tag, 64, 0, 64, REGULARITY_NODE_BUDGET,
        scope, dfs_script))
    return artifact_from_output(output)


def expected_shard_path(tile, lo_idx, scope="selected", tag=None):
    tag = tile_id(tile, scope) if tag is None else tag
    return os.path.join(HERE, "minimum_dfs_%s_p%03d.json" %
                        (tag, lo_idx))


def shard_matches(path, tile, lo_idx, hi_idx, scope="selected", n_parts=64):
    doc = load_json(path)
    if doc is None:
        return False
    spec = doc.get("spec", {})
    try:
        bounds_match = key_of_bounds(*(spec.get(name) for name in
                                       ("b0", "b1", "y0", "y1"))) \
            == tile_key(tile)
    except Exception:
        return False
    return (bounds_match and spec.get("coord") == "relative"
            and same_strip(spec) and spec.get("n_parts") == n_parts
            and spec_scope(spec) == scope
            and doc.get("lo_idx") == lo_idx
            and doc.get("hi_idx") == hi_idx)


def usable_shards(tile, scope="selected", n_parts=64, tag=None):
    """Shards on disk that this tile's merge would accept, as [lo,hi) pairs.

    A shard qualifies only if its whole spec matches the tile and requested
    root partition, and it carries no undecided leaf.  Anything else is
    treated as absent so that it is regenerated.
    """
    tag = tile_id(tile, scope) if tag is None else tag
    out = []
    for path in glob.glob(os.path.join(HERE, "minimum_dfs_%s_p*.json" % tag)):
        doc = load_json(path)
        if doc is None:
            continue
        spec = doc.get("spec", {})
        try:
            bounds = key_of_bounds(*(spec.get(name) for name in
                                     ("b0", "b1", "y0", "y1")))
        except Exception:
            continue
        if bounds != tile_key(tile) or spec.get("coord") != "relative" \
                or not same_strip(spec) or spec.get("n_parts") != n_parts \
                or spec_scope(spec) != scope:
            continue
        stats = doc.get("stats", {})
        if not isinstance(stats, dict) \
                or stats.get("counts", {}).get("UNDECIDED"):
            continue
        lo_idx, hi_idx = doc.get("lo_idx"), doc.get("hi_idx")
        if isinstance(lo_idx, int) and isinstance(hi_idx, int) \
                and 0 <= lo_idx < hi_idx <= n_parts:
            out.append((lo_idx, hi_idx))
    out.sort()
    return out


def missing_ranges(tile, scope="selected", n_parts=64, tag=None):
    """Maximal root-cell ranges that no matching usable shard covers.

    Overlapping shards would break the merge, so a shard that starts before
    the first uncovered cell is skipped rather than trusted.
    """
    out = []
    reach = 0
    for lo_idx, hi_idx in usable_shards(tile, scope, n_parts, tag):
        if lo_idx > reach:
            out.append((reach, lo_idx))
            reach = hi_idx
        elif hi_idx > reach and lo_idx == reach:
            reach = hi_idx
    if reach < n_parts:
        out.append((reach, n_parts))
    return out


def generate_regularity_sharded(tile, shard_size=SHARD_SIZE,
                                shard_budget=SHARD_NODE_BUDGET,
                                scope="selected",
                                dfs_script="minimum_dfs.py"):
    """Fill the root-cell ranges still missing, then merge in checker order.

    The budget is per INVOCATION, so ``shard_size`` is really the number of
    root cells made to share one budget.  Grouping is therefore a liability
    wherever cost is concentrated in a few cells: root cells [32,64) carry
    97.4 % of all nodes, so one pathological cell exhausts the budget and the
    invocation writes nothing, discarding every cell it had already finished.
    Generating what is missing rather than a fixed grouping means a completed
    root cell is banked permanently, whatever grouping a later run requests.
    """
    if shard_size < 1 or shard_size > 64:
        raise ValueError("shard size must be between 1 and 64 root cells")
    tag = tile_id(tile, scope)
    for gap_lo, gap_hi in missing_ranges(tile, scope):
        for lo_idx in range(gap_lo, gap_hi, shard_size):
            hi_idx = min(lo_idx + shard_size, gap_hi)
            output = run_checked(regularity_command(
                tile, tag, 64, lo_idx, hi_idx, shard_budget,
                scope, dfs_script))
            path = artifact_from_output(output)
            if not shard_matches(path, tile, lo_idx, hi_idx, scope):
                raise RuntimeError("generator emitted a mismatched shard %s" %
                                   os.path.basename(path))
    output = run_checked(command("minimum_dfs_merge.py", tag))
    merged = None
    for line in reversed(output.splitlines()):
        if " -> " in line:
            merged = os.path.join(HERE, line.split(" -> ", 1)[1].split()[0])
            break
    if merged is None:
        merged = os.path.join(HERE, "minimum_dfs_%s_merged.json" % tag)
    return merged


def generate_boundary(tile, face, scope="selected"):
    tag = tile_id(tile, scope)
    output = run_checked(command(
        "minimum_boundary_dfs.py", tile["b0"], tile["b1"],
        tile["y0"], tile["y1"], STRIP["delta"], STRIP["dmax"],
        face, tag, 64, 0, 64, 500000, "1e-6", scope))
    return artifact_from_output(output)


def replay_regularity(path, dfs_checker="minimum_dfs_check.py"):
    run_checked(command(validate_dfs_checker(dfs_checker), path))


def replay_boundary(path):
    run_checked(command("minimum_boundary_check.py", path))


def validate_tile(tile, found_regularity, found_boundaries,
                  sharded=False, shard_size=SHARD_SIZE,
                  shard_budget=SHARD_NODE_BUDGET, scope="selected",
                  dfs_script="minimum_dfs.py",
                  dfs_checker="minimum_dfs_check.py"):
    started = time.time()
    key = tile_key(tile)
    regularity = found_regularity.get(key)
    if regularity is None:
        regularity = (generate_regularity_sharded(
            tile, shard_size=shard_size, shard_budget=shard_budget,
            scope=scope, dfs_script=dfs_script)
            if sharded else generate_regularity(tile, scope, dfs_script))
    replay_regularity(regularity, dfs_checker)
    boundary_paths = {}
    for face in FACES:
        path = found_boundaries.get((key, face))
        if path is None:
            path = generate_boundary(tile, face, scope)
        replay_boundary(path)
        boundary_paths[face] = os.path.basename(path)
    return {
        "id": tile_id(tile, scope),
        "b0": float(tile["b0"]), "b1": float(tile["b1"]),
        "y0": float(tile["y0"]), "y1": float(tile["y1"]),
        "regularity": os.path.basename(regularity),
        "boundaries": boundary_paths,
        "seconds": time.time() - started,
    }


def adopt_and_validate_tile(pid, tile, sharded=False,
                            shard_size=SHARD_SIZE,
                            shard_budget=SHARD_NODE_BUDGET,
                            scope="selected",
                            dfs_script="minimum_dfs.py",
                            dfs_checker="minimum_dfs_check.py"):
    """Wait for an existing generator, then replay its complete tile."""
    with PRINT_LOCK:
        print("ADOPTED pid %d for %s" % (pid, tile_id(tile)), flush=True)
    while process_is_alive(pid):
        time.sleep(10)
    regularity, boundaries = discover(scope)
    return validate_tile(tile, regularity, boundaries, sharded,
                         shard_size, shard_budget, scope, dfs_script,
                         dfs_checker)


def state_key(tile):
    return ",".join(str(tile[name]) for name in ("b0", "b1", "y0", "y1"))


def load_state():
    doc = load_json(STATE_PATH)
    if doc is None or doc.get("format") != "minimum-atlas-sweep-v1":
        return {"format": "minimum-atlas-sweep-v1", "strip": STRIP,
                "passed": {}, "failed": {}}
    if doc.get("strip") != STRIP:
        raise RuntimeError("sweep state has a different strip")
    return doc


def write_state(state):
    temporary = STATE_PATH + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(state, handle, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(temporary, STATE_PATH)


def selected_tiles(y_min=None, lower_half=False, beta_min=None):
    return [tile for tile in tiles()
            if (not lower_half or tile["y1"] <= Decimal("-0.45"))
            and (y_min is None or tile["y0"] >= Decimal(y_min))
            and (beta_min is None or tile["b0"] >= Decimal(beta_min))]


def write_atlas(state, path, scope="selected", y_min=None, lower_half=False,
                beta_min=None):
    all_tiles = selected_tiles(y_min, lower_half, beta_min)
    if not all_tiles:
        raise RuntimeError("cannot assemble an empty atlas")
    expected = {state_key(tile) for tile in all_tiles}
    passed = set(state["passed"])
    missing = expected - passed
    if missing:
        raise RuntimeError("cannot assemble full atlas: %d tiles are not "
                           "replay-journaled" % len(missing))
    boxes = []
    for tile in all_tiles:
        item = state["passed"][state_key(tile)]
        boxes.append({
            "id": tile_id(tile, scope),
            "b0": float(tile["b0"]), "b1": float(tile["b1"]),
            "y0": float(tile["y0"]), "y1": float(tile["y1"]),
            "regularity": item["regularity"],
            "boundaries": item["boundaries"],
        })
    target = {
        "b0": float(min(tile["b0"] for tile in all_tiles)),
        "b1": float(max(tile["b1"] for tile in all_tiles)),
        "y0": float(min(tile["y0"] for tile in all_tiles)),
        "y1": float(max(tile["y1"] for tile in all_tiles)),
    }
    manifest = {
        "format": "minimum-atlas-v1",
        "scope": scope,
        "strip": STRIP,
        "target": target,
        "witness": ("minimum_witness_f4probe.json" if scope == "selected"
                    else "minimum_witness_f4probe_census.json"),
        "base_box": tile_id({"b0": Decimal("0.7"),
                             "b1": Decimal("0.8"),
                             "y0": Decimal("-0.45"),
                             "y1": Decimal("-0.35")}, scope),
        "boxes": boxes,
    }
    if scope == "selected":
        manifest["selected_count"] = 1
    else:
        manifest["zero_count"] = 2
    destination = path if os.path.isabs(path) else os.path.join(HERE, path)
    temporary = destination + ".tmp"
    with open(temporary, "w") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")
    os.replace(temporary, destination)
    print("wrote %s with %d replay-journaled boxes" %
          (os.path.basename(destination), len(boxes)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int,
                        default=max(1, (os.cpu_count() or 2) - 2))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--list", action="store_true")
    parser.add_argument("--retry-passed", action="store_true")
    parser.add_argument(
        "--sharded", action="store_true",
        help="generate regularity as independently budgeted root-cell shards")
    parser.add_argument("--shard-size", type=int, default=SHARD_SIZE)
    parser.add_argument("--shard-budget", type=int, default=SHARD_NODE_BUDGET)
    parser.add_argument(
        "--dfs-script", type=parse_dfs_script, default="minimum_dfs.py",
        help="regularity generator basename; boundaries and all replay "
             "checkers are selected independently")
    parser.add_argument(
        "--dfs-checker", type=parse_dfs_checker,
        default="minimum_dfs_check.py",
        help="regularity replay checker basename; the cached checker calls "
             "the unchanged ordinary checker after installing the exact "
             "atom cache")
    parser.add_argument(
        "--lower-half", action="store_true",
        help="schedule only the 105 cells with y <= -0.45; the other rows "
             "require the separate reflection theorem")
    parser.add_argument(
        "--adopt-running", action="store_true",
        help="reserve workers for live allowed in-directory DFS generators")
    parser.add_argument(
        "--y-min", type=str,
        help="schedule only tiles with y0 >= this.  The census strip is a "
             "WINDOW in D: where the saddle gap reaches dmax the zero leaves "
             "through gap_hi, hint fails, and the tile is unclosable at this "
             "dmax however long it is swept.  That locus is the antipodal "
             "stratum (D -> pi at y = -1/2) cut by {D = dmax}, a mirror pair "
             "about y = -1/2 sitting at y = -0.415 (beta 0.3) to -0.499 "
             "(beta 2.1).  Excluding it up front avoids paying a full "
             "regularity sweep per doomed tile, since validate_tile does "
             "regularity before boundaries.")
    parser.add_argument(
        "--scope", default="selected", choices=("selected", "census"),
        help="'census' controls EVERY angle-map zero, which is what the "
             "face-count schema's fold-freeness needs; 'selected' controls "
             "only the zeros that can be local minima")
    parser.add_argument(
        "--beta-min", type=str,
        help="schedule only tiles with b0 >= this.  For the census target "
             "with y-min=-0.45 this must be at least 0.5: the D=dmax "
             "exit curve crosses the bottom row for 0.3 <= beta < 0.5, "
             "and its location below beta=0.3 is unmeasured.")
    parser.add_argument("--write-atlas")
    args = parser.parse_args()
    # Each scope keeps its own journal: a pass under one is not a pass under
    # the other, and the tiles() order is a CONNECTED Manhattan growth from the
    # witness box, which is exactly the preconnected region the closure needs.
    global STATE_PATH
    if args.scope != "selected":
        STATE_PATH = os.path.join(
            HERE, "minimum_atlas_f4_sweep_%s_state.json" % args.scope)
    all_tiles = tiles()
    if args.list:
        print("tiles:", len(all_tiles))
        for tile in all_tiles:
            print(tile_id(tile), state_key(tile))
        return 0
    if args.write_atlas is not None:
        try:
            write_atlas(load_state(), args.write_atlas, args.scope,
                        args.y_min, args.lower_half, args.beta_min)
        except RuntimeError as exc:
            print("ATLAS NOT WRITTEN:", exc)
            return 1
        return 0
    if args.workers < 1 or (args.limit is not None and args.limit < 1):
        parser.error("workers and limit must be positive")
    if args.shard_size < 1 or 64 % args.shard_size \
            or args.shard_budget < 1:
        parser.error("shard-size must divide 64 and shard-budget is positive")
    state = load_state()
    adopted = discover_running_generators(all_tiles) \
        if args.adopt_running else {}
    if not args.retry_passed:
        adopted = {key: value for key, value in adopted.items()
                   if state_key(value[1]) not in state["passed"]}
    scheduled_tiles = selected_tiles(args.y_min, args.lower_half,
                                     args.beta_min)
    scheduled_keys = {tile_key(tile) for tile in scheduled_tiles}
    adopted = {key: value for key, value in adopted.items()
               if key in scheduled_keys}
    work = [tile for tile in scheduled_tiles
            if tile_key(tile) not in adopted
            and (args.retry_passed
                 or state_key(tile) not in state["passed"])]
    if args.limit is not None:
        work = work[:args.limit]
    regularity, boundaries = discover(args.scope)
    print("F4 sweep: %d total, %d journaled, %d adopted, %d scheduled, "
          "%d workers, generator %s, checker %s" %
          (len(scheduled_tiles),
           sum(state_key(tile) in state["passed"] for tile in scheduled_tiles),
           len(adopted), len(work),
           args.workers, args.dfs_script, args.dfs_checker), flush=True)
    failures = 0
    with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.workers) as pool:
        pending = {
            pool.submit(adopt_and_validate_tile, pid, tile,
                        args.sharded, args.shard_size,
                        args.shard_budget, args.scope,
                        args.dfs_script, args.dfs_checker): tile
            for pid, tile in adopted.values()
        }
        pending.update({
            pool.submit(validate_tile, tile, regularity, boundaries,
                        args.sharded, args.shard_size,
                        args.shard_budget, args.scope,
                        args.dfs_script, args.dfs_checker): tile
            for tile in work
        })
        for future in concurrent.futures.as_completed(pending):
            tile = pending[future]
            key = state_key(tile)
            try:
                result = future.result()
            except Exception as exc:
                failures += 1
                state["failed"][key] = str(exc)
                status = "FAILED %s: %s" % (tile_id(tile), exc)
            else:
                state["passed"][key] = result
                state["failed"].pop(key, None)
                status = "PASSED %s in %.0fs (%d/%d)" % (
                    tile_id(tile), result["seconds"],
                    sum(state_key(item) in state["passed"]
                        for item in scheduled_tiles), len(scheduled_tiles))
            with PRINT_LOCK:
                write_state(state)
                print(status, flush=True)
    print("F4 sweep finished: %d passed, %d failed" %
          (len(state["passed"]), failures), flush=True)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
