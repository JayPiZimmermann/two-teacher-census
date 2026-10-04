"""REPLAY CHECKER for a continuation certificate (`cont_cert.py`).

Same discipline as `ivcert_check.py` and `foldfree_check.py`, and for the same
reason: the search that produced the tree is prospecting and appears nowhere in
the file.  Everything except the verdict bits is REGENERATED from the spec --
the two parameter grids, the comb of walks, the parameter bisection rule, the
Newton predictor, the Krawczyk radius ladder -- and every leaf verdict is
RE-DERIVED from scratch and compared against the claim.  An undecided leaf
fails the replay.

WHAT A VALID REPLAY ESTABLISHES.  At every teacher of the regenerated grid,
for each tracked family: EXACTLY ONE zero of the census angle map lies in the
re-derived Krawczyk enclosure, the enclosures of distinct families are
pairwise disjoint as unordered student pairs, the angle-map Jacobian
determinant `separatedAngleJacDetJ` is nonzero at that zero (a sign-definite
interval enclosure over the box), and the interval Schur type is the claimed
one.  This is for the CONCRETE census map at the teachers the spec names, and
for no other teacher.

WHAT IT DOES NOT ESTABLISH -- read this before quoting the exit code.  It is
NOT the `hfold` hypothesis of `censusAngleMapJ_zeroCount_eq_witness`, which
quantifies over ALL zeros of the map at ALL teachers of the face.  A
continuation certifies fold-freeness only ALONG THE SHEETS IT TRACKS: it
cannot exclude a zero that is not a continuation of a seeded family, nor the
birth of a new pair between two grid teachers.  The bootstrap that would close
that gap -- if fold-freeness held everywhere the count would be constant, so
every zero would be a continuation -- is CIRCULAR until degenerate zeros are
excluded independently, and nothing here excludes them.  Nor does a step say
anything about the parameters between its endpoints.

Exit 0 = every leaf re-derived and matching, none undecided, families disjoint.
Run: `python3 cont_check.py cont_<tag>.json`
"""
import base64
import json
import os
import sys
import time

import cont_cert as C
import noncentered as J


def main():
    path = sys.argv[1]
    doc = json.load(open(path if os.path.isabs(path)
                         else os.path.join(C.HERE, path)))
    spec = doc["spec"]
    J.set_prec(spec.get("prec", 160))

    raw = base64.b64decode(doc["bits_b64"])
    bits = []
    for byte in raw:
        for k in range(7, -1, -1):
            bits.append((byte >> k) & 1)
    bits = bits[:doc["n_nodes"]]

    src = C.Bits(bits)
    stats = C.new_stats()
    per_family = []
    t0 = time.time()
    for seed in spec["seeds"]:
        d = C.run_family(spec, seed, src, False, stats)
        if d is None:
            break
        per_family.append(d)
    secs = time.time() - t0

    bg, yg = C.grid(spec, "b"), C.grid(spec, "y")
    print("certificate : %s" % os.path.basename(path))
    print("rectangle   : beta=[%g,%g] step %s ; y=[%g,%g] step %s "
          "(regenerated from the spec)"
          % (spec["b0"], spec["b1"],
             mp_str((bg[1] - bg[0]) if len(bg) > 1 else 0),
             spec["y0"], spec["y1"],
             mp_str((yg[1] - yg[0]) if len(yg) > 1 else 0)))
    print("grid        : %d x %d = %d teachers, %d tracked families"
          % (len(bg), len(yg), len(bg) * len(yg), len(spec["seeds"])))
    print("nodes       : %d claimed, %d consumed" % (doc["n_nodes"], src.pos))
    print("leaves      : %d re-derived and matching, %d FAILED, %d undecided"
          % (stats["ok"], stats["bad"], stats["undec"]))
    print("splits      : %d parameter bisections replayed" % stats["splits"])
    print("verdicts    : %s" % stats["types"])
    for w in stats["bad_where"][:20]:
        print("   FAILED leaf %s" % (w,))
    for w in stats["undec_where"][:20]:
        print("   UNDECIDED leaf %s" % (w,))

    if stats["bad"]:
        print("CERTIFICATE INVALID")
        return 1
    if stats["undec"]:
        print("CERTIFICATE INCOMPLETE: undecided leaves, so the continuation "
              "is not certified across the rectangle.")
        return 1
    if len(per_family) != len(spec["seeds"]):
        print("CERTIFICATE INCOMPLETE: %d of %d families completed"
              % (len(per_family), len(spec["seeds"])))
        return 1
    if src.pos != doc["n_nodes"]:
        print("CERTIFICATE INVALID: %d bits claimed, %d consumed"
              % (doc["n_nodes"], src.pos))
        return 1

    ok, clash, worst = C.disjointness(spec, per_family)
    if not ok:
        print("CERTIFICATE INVALID: families %d and %d overlap at grid "
              "teacher %s" % (clash[2], clash[3], (clash[0], clash[1])))
        return 1
    print("disjoint    : YES, worst separation %.4g (families %d,%d at grid "
          "teacher %s)" % (worst[0], worst[3], worst[4], (worst[1], worst[2])))

    # the domain bookkeeping, reported and never asserted: the fold-freeness
    # statement is about zeros of the map and does not need it, while the
    # FAMILY reading does.  The quantity is the UNORDERED gap
    # `min(D mod 2pi, 2pi - (D mod 2pi))`, not the chart's `D`: the walk is
    # free to carry the labelled representative across `D = pi` (it does, on
    # F4, because the teacher-swap mirror `y -> -1-y` makes the family exactly
    # antipodal at `y = -1/2`), and that is a change of label, not of family.
    tp = 2.0 * float(J.hi(J.PI_IV()))
    for f, d in enumerate(per_family):
        gaps = []
        for v in d.values():
            g = float(J.mid(v[1])) % tp
            gaps.append(min(g, tp - g))
        glo, ghi = min(gaps), max(gaps)
        print("family %d    : unordered gap in [%.6f, %.6f] over the grid; "
              "%s the certificate's own gap floor delta = %g"
              % (f, glo, ghi, "above" if glo > spec["delta"] else "BELOW",
                 spec["delta"]))

    print("replay      : %.1f s (%.1f ms per re-derived step)"
          % (secs, 1000.0 * secs / max(stats["ok"], 1)))
    print("VALID AND COMPLETE: at every one of the %d grid teachers, each of "
          "the %d tracked families has exactly one zero of the census angle "
          "map in its certified enclosure, the enclosures are pairwise "
          "disjoint, and `separatedAngleJacDetJ` is nonzero there."
          % (len(bg) * len(yg), len(spec["seeds"])))
    print("NOT the schema's `hfold`: fold-freeness is certified ALONG THE "
          "TRACKED SHEETS at the grid teachers only.  Zeros that are not "
          "continuations of a seeded family, and births between grid "
          "teachers, are NOT excluded here.")
    return 0


def mp_str(x):
    return "%.6g" % float(x)


if __name__ == "__main__":
    raise SystemExit(main())
