"""ROW 7a's EXCLUSION, as a replayable interval certificate on a NAMED band.

WHAT IS BEING CERTIFIED, and against which Lean statement.  The Lean theorem
`CountConstancyJ/SchurExclusionJ.censusAngleMapJ_foldFreeOn_censusStripJ` says:
on `censusStripJ seam delta dmax` with `0 < delta` and `dmax < π`, the
counting schema's `hfold` follows from three inputs —

  `CensusSchurExclusionOffFitOn S K`   no `(p, x)` with the census angle map
                                       zero and `generalJKernelTeacherSchurDetT`
                                       zero, away from the exact fit;
  `CensusTorqueAliveOffFitOn S K`      both teacher torques nonzero at every
                                       such zero;
  `CensusFoldFreeAtExactFitOn S K`     the mass-carrying input AT the exact fit.

This script discharges the first two for the CONCRETE census map on the named
rectangle, and the spec is chosen so the third is free: with `delta > b1` the
gap window starts above every teacher gap in the rectangle, so the exact fit
(gap `β`) is not in the region at all
(`censusFoldFreeAtExactFitOn_censusStripJ_of_lt_delta`).

Per box, TWO one-sided alternatives discharge the obligation, and neither
encloses a moving solution nor contracts onto one:

  MAP   some component of the angle map has an enclosure missing zero over the
        whole box, so the map cannot vanish anywhere in it.  (The row form
        `E_i = massDet(D)·G_i` is used: `E_i ≠ 0` implies `G_i ≠ 0`, which is
        the direction needed.)
  SAFE  the mean-value enclosure of `generalJKernelTeacherSchurDetT` misses
        zero AND both teacher-torque enclosures miss zero, so if the map DOES
        vanish somewhere in the box the conclusion holds there.

SCOPE, stated before any number.  This is a statement about the concrete census
map on ONE teacher rectangle and ONE gap band, never a general theorem, and the
band is chosen where the determinant test can fire: at small gaps it cannot
(measured in the README), so a band reaching down to the census gap floor is
NOT covered by this or any determinant-based certificate.

FORMAT.  The replayable bisection tree of the census certificates: the root box
is regenerated from the spec alone, the split rule is deterministic (bisect the
widest side, child 0 = lower half), the tree is walked in preorder, `1` =
internal, `0` = leaf followed by two verdict bits (MAP 0, SAFE 1, UNDECIDED 3).
An undecided leaf fails the replay.  `row7a_check.py` replays it.

Usage: python3 row7a_excl.py <b0> <b1> <y0> <y1> <delta> <dmax> <tag> [budget]
"""
import base64
import json
import os
import sys
import time

import mpmath

import census_cert as X
import face_bb_signlaw as F
import noncentered as J
import kernel_schur_mv as MV

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = {"MAP": 0, "SAFE": 1, "UNDECIDED": 3}
NAME = {v: k for k, v in CODE.items()}


def mk(lo, hi):
    return X.iv.mpf([mpmath.mpf(repr(lo)), mpmath.mpf(repr(hi))])


def root_box(spec):
    """The whole region `(beta, y, s, D)`, regenerated from the spec ALONE.

    `s = θ₁` runs one full period from the seam and `D = θ₀ − θ₁` runs over
    `[delta, dmax]` — exactly `censusStripJ seam delta dmax` in the Lean
    statement's coordinates, crossed with the teacher rectangle."""
    two_pi = float(J.hi(2 * X.PIv))
    return (spec["b0"], spec["b1"], spec["y0"], spec["y1"],
            spec["seam"], spec["seam"] + two_pi, spec["delta"], spec["dmax"])


def spec_admissible(spec):
    """The three conditions the Lean statement needs of the SPEC itself.

    `0 < delta`, `dmax < π` (the gap window stops below the antipodal locus —
    sharp, see `exists_opposite_mem_censusStripJ`), and `b1 < delta` (the gap
    window starts above every teacher gap, so the exact fit is outside the
    region and row 7b's input is discharged outright)."""
    pi_lo = float(J.lo(X.PIv))
    return (spec["delta"] > 0 and spec["dmax"] < pi_lo
            and spec["b1"] < spec["delta"] and spec["b0"] > 0)


def verdict(spec, bl, bh, yl, yh, sl, sh, dl, dh):
    B = mk(bl, bh)
    T1 = mk(sl, sh)
    Dm = mk(dl, dh)
    T0 = T1 + Dm
    Y = mk(yl, yh)
    s0, s1 = X.masses_at(Y)
    Tc = X.Teacher(B, s0, s1)
    # ORDER IS BY COST: the MAP test is the cheap one and kills almost
    # everything at coarse scales, so it goes first.  Either verdict discharges
    # the box, and the checker calls THIS function, so the orders agree by
    # construction.
    try:
        E0, E1 = X.sep_res_rowform(Tc, T0, T1, None, Dm=Dm)
        if X.sgn(E0) != 0 or X.sgn(E1) != 0:
            return "MAP"
        if any(F.rowform_centered_sign(Tc, T1, Dm, None)):
            return "MAP"
    except Exception:
        pass
    try:
        mv = MV.detT_meanvalue(B, mpmath.mpf(repr(sl)), mpmath.mpf(repr(sh)),
                               mpmath.mpf(repr(dl)), mpmath.mpf(repr(dh)))
        if X.sgn(mv) == 0:
            return None
        # the torque side conditions, on the SAME box: both students' teacher
        # torques must be enclosed away from zero for the mass-free route to
        # apply at any zero the box may hold.
        A0 = Tc.load(T0)[1]
        A1 = Tc.load(T1)[1]
        if X.sgn(A0) != 0 and X.sgn(A1) != 0:
            return "SAFE"
    except Exception:
        pass
    return None


def child(bx, bit):
    """Deterministic split: widest side, child 0 the lower half."""
    w = (bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4], bx[7] - bx[6])
    i = w.index(max(w))
    lo, hi = bx[2 * i], bx[2 * i + 1]
    m = (lo + hi) / 2
    out = list(bx)
    if bit == 0:
        out[2 * i + 1] = m
    else:
        out[2 * i] = m
    return tuple(out)


def main():
    b0, b1, y0, y1, delta, dmax = (float(sys.argv[i]) for i in range(1, 7))
    tag = sys.argv[7]
    budget = int(sys.argv[8]) if len(sys.argv) > 8 else 4000000
    spec = {"b0": b0, "b1": b1, "y0": y0, "y1": y1, "seam": 0.137,
            "delta": delta, "dmax": dmax, "minw": 1e-9, "order": "map_first"}
    if not spec_admissible(spec):
        print("SPEC REJECTED: needs 0 < b0, b1 < delta, dmax < pi")
        return 2
    print("row 7a exclusion  beta=[%g,%g] y=[%g,%g] gap=[%g,%g]  budget %d"
          % (b0, b1, y0, y1, delta, dmax, budget), flush=True)
    bits = []
    counts = {k: 0 for k in CODE}
    t0 = time.time()
    n = 0
    live = 0.0
    done_vol = 0.0

    def vol(bx):
        return ((bx[1] - bx[0]) * (bx[3] - bx[2]) * (bx[5] - bx[4])
                * (bx[7] - bx[6]))

    root = root_box(spec)
    root_vol = vol(root)
    stack = [root]
    while stack:
        bx = stack.pop()
        n += 1
        if n > budget:
            print("BUDGET EXHAUSTED at node %d (stack %d)" % (n, len(stack)),
                  flush=True)
            return 2
        v = verdict(spec, *bx)
        if v is None:
            w = max(bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4],
                    bx[7] - bx[6])
            if w < spec["minw"]:
                bits.append(0)
                bits.extend((1, 1))
                counts["UNDECIDED"] += 1
                live += vol(bx)
                done_vol += vol(bx)
                continue
            bits.append(1)
            stack.append(child(bx, 1))
            stack.append(child(bx, 0))
            continue
        bits.append(0)
        bits.extend(((CODE[v] >> 1) & 1, CODE[v] & 1))
        counts[v] += 1
        done_vol += vol(bx)
        if n % 5000 == 0:
            print("  nodes %8d  map %7d safe %6d undec %d  stack %6d  "
                  "volume decided %.6f  %.0fs"
                  % (n, counts["MAP"], counts["SAFE"], counts["UNDECIDED"],
                     len(stack), done_vol / root_vol, time.time() - t0),
                  flush=True)
    secs = time.time() - t0
    pad = (-len(bits)) % 8
    packed = bytearray()
    for i in range(0, len(bits) + pad, 8):
        byte = 0
        for k in range(8):
            byte = (byte << 1) | (bits[i + k] if i + k < len(bits) else 0)
        packed.append(byte)
    doc = {"spec": spec, "n_nodes": len(bits),
           "bits_b64": base64.b64encode(bytes(packed)).decode()}
    out = os.path.join(HERE, "row7a_%s.json" % tag)
    with open(out, "w") as fh:
        json.dump(doc, fh, separators=(",", ":"))
    print("nodes %d  map %d safe %d undecided %d" %
          (n, counts["MAP"], counts["SAFE"], counts["UNDECIDED"]))
    print("wall clock %.1f s ; %s (%.1f KB)"
          % (secs, os.path.basename(out), len(packed) / 1024.0))
    if counts["UNDECIDED"]:
        print("INCOMPLETE: %d undecided leaves (live volume %.3e)"
              % (counts["UNDECIDED"], live))
        return 1
    print("COMPLETE: every leaf discharged")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
