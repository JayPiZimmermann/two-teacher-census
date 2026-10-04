"""Selected-zero exclusion DFS on one physical census-strip face.

The face is one of ``seam_lo``, ``seam_hi``, ``gap_lo`` or ``gap_hi``.
Coordinates are always absolute, so these are the fixed physical faces of
``censusStripJ seam delta dmax`` rather than parameter-dependent faces of the
teacher-relative chart.

A leaf is accepted only when ``minimum_dfs.verdict`` excludes an angle-map
zero or proves that every possible zero is unselected (negative first pivot
or negative determinant).  ``POS_DET`` is deliberately rejected: regularity
does not imply absence from the search boundary.

Usage:
  python3 minimum_boundary_dfs.py b0 b1 y0 y1 delta dmax face tag \
      n_parts lo_idx hi_idx [budget] [minw]
"""
import base64
import json
import math
import os
import sys
import time

import census_cert as X
import minimum_dfs as M
import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))
MAP_ONLY = {"MAP_PLAIN", "MAP_CENTERED", "MAP_AUX", "MAP_COLLAR"}
# "selected" scope only has to keep SELECTED zeros off the face, so proving a
# zero unselected is a legitimate exclusion.  "census" scope feeds `hint` of
# `censusAngleMapJ_zeroCount_eq_witness`, which demands the face carry NO zero
# of the angle map at all; a negative pivot or determinant says a zero is a
# saddle, not that it is absent, so only the map exclusions survive.
EXCLUDING = MAP_ONLY | {"NEG_T00", "NEG_DET"}


def excluding_for(spec):
    return MAP_ONLY if spec.get("scope") == "census" else EXCLUDING


def root_box(spec):
    two_pi = float(J.hi(2 * X.PIv))
    face = spec["face"]
    if face == "seam_lo":
        return (spec["b0"], spec["b1"], spec["y0"], spec["y1"],
                spec["seam"], spec["seam"], spec["delta"], spec["dmax"])
    if face == "seam_hi":
        seam_hi = spec["seam"] + two_pi
        return (spec["b0"], spec["b1"], spec["y0"], spec["y1"],
                seam_hi, seam_hi, spec["delta"], spec["dmax"])
    if face == "gap_lo":
        return (spec["b0"], spec["b1"], spec["y0"], spec["y1"],
                spec["seam"], spec["seam"] + two_pi,
                spec["delta"], spec["delta"])
    if face == "gap_hi":
        return (spec["b0"], spec["b1"], spec["y0"], spec["y1"],
                spec["seam"], spec["seam"] + two_pi,
                spec["dmax"], spec["dmax"])
    raise ValueError("unknown boundary face")


def root_cells(spec):
    root = root_box(spec)
    n = spec["n_parts"]
    axis = 3 if spec["face"].startswith("seam") else 2
    lo, hi = root[2 * axis], root[2 * axis + 1]
    cells = []
    for i in range(n):
        out = list(root)
        out[2 * axis] = lo + (hi - lo) * i / n
        out[2 * axis + 1] = lo + (hi - lo) * (i + 1) / n
        cells.append(tuple(out))
    return cells


def partition_covers(spec):
    root, cells = root_box(spec), root_cells(spec)
    axis = 3 if spec["face"].startswith("seam") else 2
    if not cells or cells[0][2 * axis] != root[2 * axis] \
            or cells[-1][2 * axis + 1] != root[2 * axis + 1]:
        return False
    for i, cell in enumerate(cells):
        for j in range(4):
            if j != axis and cell[2 * j:2 * j + 2] != root[2 * j:2 * j + 2]:
                return False
        if i and cells[i - 1][2 * axis + 1] != cell[2 * axis]:
            return False
    return True


def child(spec, box, bit):
    return M.child(spec, box, bit)


def excludes(spec, box):
    return M.verdict(spec, *box) in excluding_for(spec)


def encode(bits):
    padding = (-len(bits)) % 8
    packed = bytearray()
    for i in range(0, len(bits) + padding, 8):
        byte = 0
        for k in range(8):
            byte = (byte << 1) | (bits[i + k]
                                  if i + k < len(bits) else 0)
        packed.append(byte)
    return base64.b64encode(bytes(packed)).decode()


def run(spec, lo_idx, hi_idx, budget):
    roots = root_cells(spec)
    bits = []
    nodes = leaves = 0
    started = time.time()
    for idx in range(lo_idx, hi_idx):
        stack = [roots[idx]]
        while stack:
            box = stack.pop()
            nodes += 1
            if nodes > budget:
                raise RuntimeError("budget exhausted at cell %d" % idx)
            if excludes(spec, box):
                bits.append(0)
                leaves += 1
                continue
            width = max(box[2 * i + 1] - box[2 * i] for i in range(4))
            if width < spec["minw"]:
                raise RuntimeError("unresolved boundary box %s" % (box,))
            bits.append(1)
            stack.append(child(spec, box, 1))
            stack.append(child(spec, box, 0))
            if nodes % 5000 == 0:
                print("  face %s cell %d nodes %d leaves %d stack %d %.0fs" %
                      (spec["face"], idx, nodes, leaves, len(stack),
                       time.time() - started), flush=True)
    return {"spec": spec, "lo_idx": lo_idx, "hi_idx": hi_idx,
            "n_nodes": len(bits), "bits_b64": encode(bits),
            "stats": {"visited_nodes": nodes, "leaves": leaves,
                      "seconds": time.time() - started}}


def main():
    if len(sys.argv) < 12:
        print(__doc__)
        return 2
    b0, b1, y0, y1, delta, dmax = map(float, sys.argv[1:7])
    face, tag = sys.argv[7:9]
    n_parts, lo_idx, hi_idx = map(int, sys.argv[9:12])
    budget = int(sys.argv[12]) if len(sys.argv) > 12 else 2_000_000
    minw = float(sys.argv[13]) if len(sys.argv) > 13 else 1e-8
    scope = sys.argv[14] if len(sys.argv) > 14 else "selected"
    spec = {"b0": b0, "b1": b1, "y0": y0, "y1": y1,
            "seam": 0.137, "delta": delta, "dmax": dmax,
            "face": face, "n_parts": n_parts, "minw": minw,
            "coord": "absolute", "row_mv_maxw": 0.1,
            "kraw_maxw": 0.1, "mv_maxw": 0.05,
            "split": "widest"}
    if scope != "selected":
        spec["scope"] = scope
    if not (0 < b0 <= b1 < math.pi and -1 < y0 <= y1 < 1
            and 0 < delta < dmax < math.pi
            and face in {"seam_lo", "seam_hi", "gap_lo", "gap_hi"}
            and scope in ("selected", "census")
            and 0 <= lo_idx < hi_idx <= n_parts
            and partition_covers(spec)):
        print("SPEC REJECTED")
        return 2
    try:
        doc = run(spec, lo_idx, hi_idx, budget)
    except RuntimeError as exc:
        print(exc)
        return 1
    path = os.path.join(HERE, "minimum_boundary_%s_%s_p%03d.json" %
                        (tag, face, lo_idx))
    with open(path, "w") as handle:
        json.dump(doc, handle, separators=(",", ":"))
    print(json.dumps(doc["stats"], indent=1))
    print("wrote", os.path.basename(path))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
