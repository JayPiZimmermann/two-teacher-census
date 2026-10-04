"""Validate and optionally replay a connected minimum/census box atlas.

Each atlas box names one full-period relative regularity forest and the four
absolute physical-boundary forests required by Lean's box theorem. A base
witness manifest supplies the claimed selected count. Closed boxes that
overlap form a chain used by the full replay to transport that count.
Only the full replay checks the count claim; ``--structure-only`` checks the
manifest graph and coverage without establishing the claimed count.
``CERTIFICATE.md`` (status 2026-08-24) is the authoritative claim ledger.

Default mode independently replays every named artifact.  ``--structure-only``
checks specs, full-period coverage, the four boundary names, and connectivity
without repeating the potentially long interval arithmetic.  ``--workers``
parallelizes independent artifact replays without changing their semantics.

Usage: python3 minimum_atlas_check.py minimum_atlas_<tag>.json
       python3 minimum_atlas_check.py minimum_atlas_<tag>.json --structure-only
       python3 minimum_atlas_check.py minimum_atlas_<tag>.json --workers 20
       python3 minimum_atlas_check.py minimum_atlas_<tag>.json --workers 20 \
           --dfs-checker minimum_dfs_check_cached.py
       MINIMUM_REPLAY_WORKERS=4 python3 minimum_atlas_check.py \
           minimum_atlas_<tag>.json --workers 4 \
           --dfs-checker minimum_dfs_check_parallel.py
       MINIMUM_REPLAY_WORKERS=4 python3 minimum_atlas_check.py \
           minimum_atlas_<tag>.json --workers 4 \
           --dfs-checker minimum_dfs_check_dynamic.py
"""
import argparse
import concurrent.futures
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FACES = {"seam_lo", "seam_hi", "gap_lo", "gap_hi"}
DFS_CHECKERS = ("minimum_dfs_check.py", "minimum_dfs_check_cached.py",
                "minimum_dfs_check_parallel.py",
                "minimum_dfs_check_dynamic.py")


def local_path(path):
    return path if os.path.isabs(path) else os.path.join(HERE, path)


def load(path):
    with open(local_path(path)) as handle:
        return json.load(handle)


def same_box(spec, box, strip):
    keys = ("b0", "b1", "y0", "y1")
    return (all(spec.get(k) == box[k] for k in keys)
            and spec.get("seam") == strip["seam"]
            and spec.get("delta") == strip["delta"]
            and spec.get("dmax") == strip["dmax"])


def overlaps(a, b):
    return (max(a["b0"], b["b0"]) <= min(a["b1"], b["b1"])
            and max(a["y0"], b["y0"]) <= min(a["y1"], b["y1"]))


def validate_target(boxes, target):
    """Prove by a finite endpoint grid that the closed boxes cover target."""
    for key in ("b0", "b1", "y0", "y1"):
        if key not in target:
            raise ValueError("target rectangle is missing %s" % key)
    if not (target["b0"] < target["b1"]
            and target["y0"] < target["y1"]):
        raise ValueError("target rectangle is degenerate")
    for box in boxes:
        if not (target["b0"] <= box["b0"] <= box["b1"] <= target["b1"]
                and target["y0"] <= box["y0"] <= box["y1"]
                <= target["y1"]):
            raise ValueError("box %s is not contained in target" % box["id"])
    xs = sorted({target["b0"], target["b1"]}
                | {value for box in boxes for value in (box["b0"], box["b1"])
                   if target["b0"] <= value <= target["b1"]})
    ys = sorted({target["y0"], target["y1"]}
                | {value for box in boxes for value in (box["y0"], box["y1"])
                   if target["y0"] <= value <= target["y1"]})
    for x0, x1 in zip(xs, xs[1:]):
        for y0, y1 in zip(ys, ys[1:]):
            if x1 <= target["b0"] or target["b1"] <= x0 \
                    or y1 <= target["y0"] or target["y1"] <= y0:
                continue
            if not any(box["b0"] <= x0 and x1 <= box["b1"]
                       and box["y0"] <= y0 and y1 <= box["y1"]
                       for box in boxes):
                raise ValueError("atlas misses target cell [%s,%s] x [%s,%s]"
                                 % (x0, x1, y0, y1))


def spec_scope(spec):
    """Old artifacts predate explicit scopes and are selected-only."""
    return spec.get("scope", "selected")


def validate_box(box, strip, scope):
    reg = load(box["regularity"])
    spec = reg["spec"]
    if (not same_box(spec, box, strip) or spec.get("coord") != "relative"
            or spec_scope(spec) != scope):
        raise ValueError("regularity spec mismatch for %s" % box["id"])
    if reg["lo_idx"] != 0 or reg["hi_idx"] != spec["n_parts"]:
        raise ValueError("regularity forest is sharded/incomplete for %s"
                         % box["id"])
    boundaries = box.get("boundaries", {})
    if set(boundaries) != FACES:
        raise ValueError("four physical boundary faces not named for %s"
                         % box["id"])
    for face, path in boundaries.items():
        doc = load(path)
        bspec = doc["spec"]
        if (not same_box(bspec, box, strip)
                or bspec.get("coord") != "absolute"
                or bspec.get("face") != face
                or spec_scope(bspec) != scope):
            raise ValueError("boundary spec mismatch for %s/%s"
                             % (box["id"], face))
        if doc["lo_idx"] != 0 or doc["hi_idx"] != bspec["n_parts"]:
            raise ValueError("boundary forest is incomplete for %s/%s"
                             % (box["id"], face))


def connected_order(boxes, base_id):
    by_id = {box["id"]: box for box in boxes}
    if len(by_id) != len(boxes) or base_id not in by_id:
        raise ValueError("box ids are duplicated or base_box is absent")
    seen = {base_id}
    order = [base_id]
    while len(seen) < len(boxes):
        added = None
        for box in boxes:
            if box["id"] in seen:
                continue
            if any(overlaps(box, by_id[old]) for old in seen):
                added = box["id"]
                break
        if added is None:
            raise ValueError("box overlap graph is disconnected")
        seen.add(added)
        order.append(added)
    return order


def replay(script, artifact):
    result = subprocess.run(
        [sys.executable, os.path.join(HERE, script), artifact], cwd=HERE,
        text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        tail = "\n".join(result.stdout.splitlines()[-20:])
        raise ValueError("replay failed: %s %s\n%s" %
                         (script, artifact, tail))
    return result.stdout


def replay_all(jobs, workers):
    if workers == 1:
        for script, artifact in jobs:
            print(replay(script, artifact), end="")
        return
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        pending = {pool.submit(replay, script, artifact): (script, artifact)
                   for script, artifact in jobs}
        for future in concurrent.futures.as_completed(pending):
            script, artifact = pending[future]
            output = future.result()
            print("--- replayed %s %s ---" % (script, artifact))
            print(output, end="")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest")
    parser.add_argument("--structure-only", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument(
        "--dfs-checker", choices=DFS_CHECKERS,
        default="minimum_dfs_check.py",
        help="fresh-process regularity replay entry point; the cached choice "
             "calls the same checker logic with exact atom memoization, and "
             "the parallel/dynamic choices use MINIMUM_REPLAY_WORKERS fresh "
             "workers")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    manifest = load(args.manifest)
    try:
        if manifest.get("format") != "minimum-atlas-v1":
            raise ValueError("unknown atlas format")
        strip = manifest["strip"]
        scope = manifest.get("scope", "selected")
        if scope == "selected":
            claimed_count = manifest["selected_count"]
            witness_key = "expected_selected"
            count_name = "selected count"
        elif scope == "census":
            claimed_count = manifest["zero_count"]
            witness_key = "expected_zeros"
            count_name = "angle-map zero count"
        else:
            raise ValueError("unknown atlas scope %s" % scope)
        boxes = manifest["boxes"]
        if not boxes:
            raise ValueError("atlas has no boxes")
        witness = load(manifest["witness"])
        if (witness.get("format") != "minimum-witness-v1"
                or witness.get("scope", "selected") != scope
                or witness["strip"] != strip
                or witness[witness_key] != claimed_count):
            raise ValueError("base witness does not match atlas claim")
        for box in boxes:
            validate_box(box, strip, scope)
        order = connected_order(boxes, manifest["base_box"])
        if "target" in manifest:
            validate_target(boxes, manifest["target"])
        base = next(box for box in boxes if box["id"] == manifest["base_box"])
        teacher = witness["teacher"]
        if not (base["b0"] <= teacher["beta"] <= base["b1"]
                and base["y0"] <= teacher["y"] <= base["y1"]):
            raise ValueError("witness teacher is outside the base box")
        if not args.structure_only:
            jobs = [("minimum_witness_check.py", manifest["witness"])]
            for box_id in order:
                box = next(item for item in boxes if item["id"] == box_id)
                jobs.append((args.dfs_checker, box["regularity"]))
                for face in sorted(FACES):
                    jobs.append(("minimum_boundary_check.py",
                                 box["boundaries"][face]))
            replay_all(jobs, args.workers)
    except (KeyError, ValueError, OSError) as exc:
        print("ATLAS INVALID:", exc)
        return 1
    print("atlas        : %s" % os.path.basename(args.manifest))
    print("boxes        : %d, connected order %s" %
          (len(boxes), " -> ".join(order)))
    print("strip        : seam=%s D=[%s,%s]" %
          (strip["seam"], strip["delta"], strip["dmax"]))
    print("scope        : %s" % scope)
    if "target" in manifest:
        target = manifest["target"]
        print("target       : beta=[%s,%s], y=[%s,%s], fully covered" %
              (target["b0"], target["b1"], target["y0"], target["y1"]))
    if args.structure_only:
        print("ATLAS STRUCTURE VALID (REPLAY SKIPPED): manifest claims "
              "%s = %d; the count was not independently replayed" %
              (count_name, claimed_count))
    else:
        print("ATLAS VALID AND REPLAYED: %s = %d on every named box" %
              (count_name, claimed_count))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
