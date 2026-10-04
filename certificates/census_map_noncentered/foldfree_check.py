"""REPLAY CHECKER for a fold-freeness certificate.

Same discipline as `ivcert_check.py`, and for the same reason: the search that
produced the tree is prospecting and appears nowhere in the file.  The root box
is REGENERATED from the spec, every box is reconstructed from the deterministic
split rule (bisect the widest side, child 0 = lower half), and every leaf
verdict is RE-DERIVED from scratch.  An undecided leaf fails the replay.

WHAT A VALID REPLAY ESTABLISHES.  Coverage is structural: a node is a leaf or
the union of its two children, so the leaves tile the root box.  If every leaf
is valid then at no point of the rectangle can the census angle map and the
mass-free type-boundary function vanish together -- except inside the declared
exact-fit tube, where the mass-free function is identically zero by a theorem
of the tree and the schema's OTHER route applies.  That is exactly the `hfold`
hypothesis of `censusAngleMapJ_familyCount_eq_witness_on_censusBox`, for the
CONCRETE census map on the named rectangle and for no other teacher.

Exit 0 = every leaf verified AND none undecided.
Run: `python3 foldfree_check.py foldfree_<tag>.json`
"""
import base64
import json
import os
import sys

import foldfree_sweep as S


def main():
    path = sys.argv[1]
    doc = json.load(open(path if os.path.isabs(path)
                         else os.path.join(S.HERE, path)))
    spec = doc["spec"]
    raw = base64.b64decode(doc["bits_b64"])
    bits = []
    for byte in raw:
        for k in range(7, -1, -1):
            bits.append((byte >> k) & 1)
    bits = bits[:doc["n_nodes"]]
    if not S.partition_covers(spec):
        print("SPEC REJECTED: the regenerated partition does not tile the "
              "region, so coverage is not structural and no completeness "
              "claim follows.")
        return 1
    cells = S.root_cells(spec)
    print("partition   : %d root cells regenerated from the spec, union "
          "verified to be the region exactly" % len(cells))

    pos = [0]
    stats = {"ok": 0, "bad": 0, "undec": 0}
    bad = []

    def walk(bx, depth):
        if pos[0] >= len(bits):
            bad.append(("bitstream exhausted", bx))
            stats["bad"] += 1
            return
        b = bits[pos[0]]
        pos[0] += 1
        if b == 1:
            if depth > 400:
                bad.append(("depth limit", bx))
                stats["bad"] += 1
                return
            walk(S.child(bx, 0), depth + 1)
            walk(S.child(bx, 1), depth + 1)
            return
        code = (bits[pos[0]] << 1) | bits[pos[0] + 1]
        pos[0] += 2
        claimed = S.NAME[code]
        if claimed == "UNDECIDED":
            stats["undec"] += 1
            return
        got = S.verdict(spec, *bx)
        if claimed == "FIT":
            ok = S.inside_fit_tube(spec, *bx)
        else:
            ok = (got == claimed)
        if ok:
            stats["ok"] += 1
        else:
            stats["bad"] += 1
            if len(bad) < 20:
                bad.append((claimed, bx))

    sys.setrecursionlimit(100000)
    for c in cells:
        walk(c, 0)

    print("certificate : %s" % os.path.basename(path))
    print("rectangle   : beta=[%g,%g] y=[%g,%g] (regenerated from the spec)"
          % (spec["b0"], spec["b1"], spec["y0"], spec["y1"]))
    print("nodes       : %d claimed, %d consumed" % (doc["n_nodes"], pos[0]))
    print("test order  : %s (recorded in the spec; either verdict discharges "
          "a box, so the order is load-bearing for REPRODUCIBILITY only)"
          % spec.get("order", "unspecified"))
    print("leaves      : %d verified, %d FAILED, %d undecided"
          % (stats["ok"], stats["bad"], stats["undec"]))
    for b in bad:
        print("   FAILED leaf %s %s" % b)
    if stats["bad"]:
        print("CERTIFICATE INVALID")
        return 1
    if stats["undec"]:
        print("CERTIFICATE INCOMPLETE: undecided leaves, so fold-freeness is "
              "NOT discharged on this rectangle.")
        return 1
    print("FOLD-FREENESS DISCHARGED on this rectangle (outside the declared "
          "exact-fit tube): the census angle map and the mass-free "
          "type-boundary function never vanish together here.  This is the "
          "`hfold` hypothesis of the Lean face theorem, for THIS map and this "
          "rectangle only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
