"""REPLAY CHECKER for a row-7a exclusion certificate.

Reads `row7a_<tag>.json` and re-derives everything except the bitstring: the
root box is regenerated from the SPEC alone, every child box from the
deterministic split rule, and every leaf's verdict is RECOMPUTED with
`row7a_excl.verdict`.  A claimed verdict that does not recompute fails the
replay, and so does any undecided leaf, any leftover bit, and any spec that
does not meet the Lean statement's own conditions on the region
(`0 < delta`, `dmax < π`, `b1 < delta`).

WHAT THIS IS AND IS NOT (conventions §0f).  This is a REPLAY, not an
independent verification: the checker calls the emitter's own verdict
function, so it re-establishes that the stored tree is the tree that function
produces on this spec — it does not cross-check the interval arithmetic
against a second implementation.  That is the same standing of every replay
checker in this directory, and it is stated here rather than left to be
assumed.

Coverage is STRUCTURAL, not asserted: an internal node's two children are the
two halves of it about the midpoint of its widest side, so the leaves of the
walked tree tile the root box exactly, and the root box is the region named in
the spec.

Usage: python3 row7a_check.py row7a_<tag>.json
"""
import base64
import json
import os
import sys

import noncentered as J
import census_cert as X
import row7a_excl as E


def main():
    path = sys.argv[1]
    if not os.path.isabs(path):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
    doc = json.load(open(path))
    spec = doc["spec"]
    print("certificate : %s" % os.path.basename(path))
    print("spec        : beta=[%g,%g] y=[%g,%g] gap=[%g,%g] seam=%g minw=%g"
          % (spec["b0"], spec["b1"], spec["y0"], spec["y1"], spec["delta"],
             spec["dmax"], spec["seam"], spec["minw"]))
    if not E.spec_admissible(spec):
        print("REPLAY FAILED: the spec does not meet the Lean statement's "
              "conditions (0 < b0, b1 < delta, dmax < pi)")
        return 1
    raw = base64.b64decode(doc["bits_b64"])
    bits = []
    for byte in raw:
        for k in range(7, -1, -1):
            bits.append((byte >> k) & 1)
    n_nodes = doc["n_nodes"]
    root = E.root_box(spec)
    print("root box    : regenerated from the spec, "
          "s=[%.6f,%.6f] D=[%g,%g]" % (root[4], root[5], root[6], root[7]))

    pos = 0
    leaves = 0
    failed = 0
    undecided = 0
    counts = {"MAP": 0, "SAFE": 0}
    vol_leaf = 0.0

    def vol(bx):
        return ((bx[1] - bx[0]) * (bx[3] - bx[2]) * (bx[5] - bx[4])
                * (bx[7] - bx[6]))

    stack = [root]
    while stack:
        bx = stack.pop()
        if pos >= len(bits):
            print("REPLAY FAILED: bitstream exhausted at node %d" % pos)
            return 1
        tag = bits[pos]
        pos += 1
        if tag == 1:
            stack.append(E.child(bx, 1))
            stack.append(E.child(bx, 0))
            continue
        if pos + 1 >= len(bits) + 1:
            print("REPLAY FAILED: truncated verdict bits")
            return 1
        code = (bits[pos] << 1) | bits[pos + 1]
        pos += 2
        leaves += 1
        vol_leaf += vol(bx)
        if code == 3:
            undecided += 1
            continue
        want = E.NAME.get(code)
        got = E.verdict(spec, *bx)
        if got != want:
            failed += 1
            if failed <= 5:
                print("  LEAF MISMATCH: claimed %s, recomputed %s, box %s"
                      % (want, got, bx))
        else:
            counts[want] += 1

    print("bits        : %d claimed, %d consumed (the stored count is the "
          "BIT count: one bit per internal node, three per leaf)"
          % (n_nodes, pos))
    print("leaves      : %d ; map %d ; safe %d ; FAILED %d ; undecided %d"
          % (leaves, counts["MAP"], counts["SAFE"], failed, undecided))
    print("coverage    : leaf volume %.12e against root volume %.12e"
          % (vol_leaf, vol(root)))
    if pos != n_nodes:
        print("REPLAY FAILED: node count mismatch")
        return 1
    if failed or undecided:
        print("REPLAY FAILED")
        return 1
    print("CERTIFICATE VALID AND COMPLETE: every leaf recomputes, none "
          "undecided, so on this teacher rectangle and this gap band the "
          "census angle map has no zero at which the mass-free type-boundary "
          "function vanishes, and both teacher torques are alive at every "
          "zero -- i.e. CensusSchurExclusionOffFitOn and "
          "CensusTorqueAliveOffFitOn hold there.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
