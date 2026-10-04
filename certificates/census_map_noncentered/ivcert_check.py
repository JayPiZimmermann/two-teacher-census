"""REPLAY CHECKER for a separated-enumeration completeness certificate.

Owner artifact standard (2026-08-08): inside the completeness certificate the
only evidence classes are Lean, iv certificates, and the README glue.  A
search log is not evidence.  This file is the other half of the iv class for
this lane: `face_bb_signlaw.py` SEARCHES (prospecting -- heap, priorities,
orderings, none of it evidence) and emits `ivcert_<tag>.json`; this checker
REPLAYS that certificate using nothing from the search.

WHAT THE CERTIFICATE IS.  A bisection tree over the `(th1, D)` domain.  Its
root cells are determined by the spec (`beta, y, delta, seam, dmax_pad,
atom_eps, cell`) and regenerated here, never read from the certificate.  Each
leaf is named by its PATH -- the sequence of bisection choices from its root
cell -- and carries a VERDICT.  Paths rather than coordinates, because the
split rule is deterministic (bisect the wider side, `0` = lower half), so the
checker reconstructs every box exactly and no float comparison is needed to
identify a box.

WHY A CHECKED TREE IS A COMPLETENESS PROOF.  Coverage is structural: a node is
either a leaf or the union of its two children, so the leaves of a finite tree
tile the root cells exactly, and the root cells tile the domain.  If every leaf
carries a verdict that is VALID -- meaning the leaf provably holds no solution,
or holds exactly one and is enclosed -- then the domain holds exactly the
enclosed solutions and nothing else.  That is the completeness statement, and
it is what this script verifies leaf by leaf.

THE VERDICTS AND WHAT EACH ONE CLAIMS (each re-derived here from scratch):

  PREFILTER        the sign law excludes the box.  Re-derived by the order
                   comparisons of `face_bb_signlaw.prefilter_kills`, whose
                   soundness is the Lean pair
                   `separatedPairCritical_student0/1_offsets_straddle`
                   (positive masses) and `..._offsets_same_side` (mixed).
  EXCL_PLAIN       the outward-rounded enclosure of `G0` or of `G1` over the
                   box does not contain zero, so no solution lies in it.
  EXCL_CENTERED    the same, for the mean-value enclosure
                   `G(X) subset G(xc) + J(X).(X - xc)`, in EITHER of the two
                   admissible forms: the Cramer one, or the division-free row
                   form `E = massDet * G` whose derivative does not divide by
                   massDet and is ~180x tighter.  Both enclose the same
                   quantity up to a positive factor, so accepting either is
                   sound and keeps older certificates replayable.
  KRAWCZYK_EMPTY   `K(X) cap X = empty`, so no solution in the box.
  KRAWCZYK_UNIQUE  `K(X) subset int(X)`, so EXACTLY ONE solution in the box.
  KRAWCZYK_INFLATED  the same on an inflated box `Z' ⊇ Z`; the leaf then holds
                   at most that one solution.
  UNDECIDED_*      no verdict.  A certificate containing any of these is NOT
                   a completeness proof, and the checker says so.

SCOPE (skill 1c).  This is an interval-arithmetic certificate about the
CONCRETE census map at ONE teacher.  It is never the evidence for a general
theorem, and it says nothing about any other teacher: transporting it across a
face is the count-constancy schema's job, from its own inputs.

Arithmetic: `mpmath.iv` with directed rounding at the working precision of
`census_cert`, so every enclosure printed by this checker contains the true
value.  Replay: `python3 ivcert_check.py ivcert_<tag>.json`.
Exit 0 = every leaf verified AND no undecided leaf.
"""
import json
import math
import os
import sys

import mpmath

import census_cert as X
import face_bb_signlaw as F
import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))


def root_cells(spec):
    """Regenerate the root cells from the spec ALONE (never from the file)."""
    beta_f, y_f = spec["beta"], spec["y"]
    delta = mpmath.mpf(repr(spec["delta"]))
    seam = mpmath.mpf(repr(spec["seam"]))
    pad = mpmath.mpf(repr(spec["dmax_pad"]))
    aeps = mpmath.mpf(repr(spec["atom_eps"]))
    cell = spec["cell"]
    B = (J.I(str(beta_f)) if spec.get("beta_hi") is None else
         X.iv.mpf([mpmath.mpf(str(beta_f)), mpmath.mpf(str(spec["beta_hi"]))]))
    Y = (J.I(str(y_f)) if spec.get("y_hi") is None else
         X.iv.mpf([mpmath.mpf(str(y_f)), mpmath.mpf(str(spec["y_hi"]))]))
    s0, s1 = X.masses_at(Y)
    Tc = X.Teacher(B, s0, s1)

    two_pi = J.hi(2 * X.PIv)
    dpi = J.hi(X.PIv)
    s_lo, s_hi = seam, seam + two_pi
    d_hi = dpi + pad
    a, b = J.endpoints(Tc.B)

    cuts = {s_lo, s_hi}
    for z in (a, b, a + dpi, b + dpi, a - dpi, b - dpi,
              a + two_pi, b + two_pi, dpi + seam):
        if s_lo < z < s_hi:
            cuts.add(z)
    cuts = {z for z in cuts if z in (s_lo, s_hi)} | {
        w for z in cuts if z not in (s_lo, s_hi)
        for w in (z - aeps, z + aeps) if s_lo < w < s_hi}
    cuts.add(s_lo)
    cuts.add(s_hi)

    def dedupe(cs, tol=mpmath.mpf("1e-12")):
        cs = sorted(cs)
        out = [cs[0]]
        for z in cs[1:]:
            if z - out[-1] > tol:
                out.append(z)
        return out

    def slice_up(cs):
        out = []
        for i in range(len(cs) - 1):
            u, v = cs[i], cs[i + 1]
            n = max(1, int(math.ceil(float(v - u) / cell)))
            for k in range(n):
                lo, hi = u + (v - u) * k / n, u + (v - u) * (k + 1) / n
                if hi > lo:
                    out.append((lo, hi))
        return out

    sf = slice_up(dedupe(cuts))
    dcuts = list({delta, d_hi} | {z for z in (a, b, a + dpi, b + dpi,
                                              two_pi - b, two_pi - a)
                                  if delta < z < d_hi})
    dcuts = [z for z in dcuts if z in (delta, d_hi)] + [
        w for z in dcuts if z not in (delta, d_hi)
        for w in (z - aeps, z + aeps) if delta < w < d_hi]
    df = slice_up(dedupe(dcuts))
    return Tc, [(x[0], x[1], z[0], z[1]) for x in sf for z in df]


def child(bx, bit):
    """The deterministic split: bisect the wider side, bit 0 = lower half."""
    s0_, s1_, d0_, d1_ = bx
    if s1_ - s0_ >= d1_ - d0_:
        m = (s0_ + s1_) / 2
        return (s0_, m, d0_, d1_) if bit == 0 else (m, s1_, d0_, d1_)
    m = (d0_ + d1_) / 2
    return (s0_, s1_, d0_, m) if bit == 0 else (s0_, s1_, m, d1_)


def branch_data(Tc, bx):
    T1 = X.iv.mpf([bx[0], bx[1]])
    Dm = X.iv.mpf([bx[2], bx[3]])
    T0 = T1 + Dm
    bb, ok = {}, True
    for key, Xi in (("n0", T0), ("n1", T1), ("n0b", T0 - Tc.B),
                    ("n1b", T1 - Tc.B), ("nD", Dm)):
        ps = J.branch_pieces(Xi)
        if len(ps) != 1:
            ok = False
            break
        bb[key] = ps[0][1]
    return T1, Dm, (bb if ok else None)


def verify_leaf(Tc, bx, verdict, beta_lo, beta_hi, positive):
    """Re-derive the claimed verdict from scratch.

    Returns `(ok, encl)`: `ok` says whether the claimed verdict is re-derived,
    and `encl` is the `(s, D)` enclosure of the solution for the two verdicts
    that assert one, `None` otherwise.  The enclosure is what the boundary
    containment report below consumes; it is re-derived here, never read from
    the certificate."""
    if verdict == "PREFILTER":
        return F.prefilter_kills(float(bx[0]), float(bx[1]), float(bx[2]),
                                 float(bx[3]), beta_lo, positive,
                                 beta_hi), None
    T1, Dm, br = branch_data(Tc, bx)
    if verdict == "EXCL_PLAIN":
        r = X.sep_res(Tc, T1 + Dm, T1, br, Dm=Dm)
        return (r is not None
                and (X.sgn(r[0]) != 0 or X.sgn(r[1]) != 0)), None
    if verdict == "EXCL_CENTERED":
        # TWO mean-value forms are admissible here and the leaf is valid if
        # EITHER signs a component, because both are enclosures of the same
        # quantity up to the positive factor massDet(D): the Cramer form
        # (`G_i`, the original) and the division-free row form
        # (`E_i = massDet * G_i`), whose derivative never divides by massDet
        # and is therefore ~180x tighter in the median.  Accepting both keeps
        # every certificate emitted before the row form existed replayable,
        # and lets the ones emitted after it use the sharper test.
        return (any(F.sep_res_centered_sign(Tc, T1, Dm, br))
                or any(F.rowform_centered_sign(Tc, T1, Dm, br))), None
    if verdict in ("KRAWCZYK_EMPTY", "KRAWCZYK_UNIQUE"):
        k, N0, N1 = X._krawczyk_sD(Tc, T1, Dm, br)
        if verdict == "KRAWCZYK_EMPTY":
            return k == "empty", None
        if k != "unique":
            return False, None
        # Iterate the certified contraction: every step is another Krawczyk
        # `unique` verdict on the previous enclosure, so the refined box is
        # still a proof, and the sharper enclosure is what the boundary margin
        # and the lattice test below need.
        S, DD = N0, N1
        for _ in range(60):
            kk, M0, M1 = X._krawczyk_sD(Tc, S, DD, br)
            if kk != "unique" or (J.width(M0) >= J.width(S)
                                  and J.width(M1) >= J.width(DD)):
                break
            S, DD = M0, M1
        return True, (S, DD)
    if verdict == "KRAWCZYK_INFLATED":
        z = F._certify_inflated(Tc, bx, br)
        return z is not None, z
    return False, None


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    if path is None:
        print(__doc__)
        sys.exit(2)
    doc = json.load(open(path if os.path.isabs(path)
                         else os.path.join(HERE, path)))
    spec = doc["spec"]
    if "bits_b64" in doc:
        # PACKED form: the tree arrives as a preorder bitstring.  Decoding is
        # the exact inverse of ivcert_pack.pack and uses nothing else from the
        # file -- the boxes are still regenerated from the spec below.
        import ivcert_pack as P
        leaves = P.unpack(doc["bits_b64"], doc["n_nodes"],
                          spec["n_root_cells"])
    else:
        leaves = doc["leaves"]
    Tc, roots = root_cells(spec)
    if len(roots) != spec["n_root_cells"]:
        print("ROOT CELLS DISAGREE: regenerated %d, certificate claims %d"
              % (len(roots), spec["n_root_cells"]))
        sys.exit(1)
    beta_lo = float(J.lo(Tc.B))
    beta_hi = float(J.hi(Tc.B))
    # THE SIGN-LAW BRANCH MUST COME FROM THE WHOLE PARAMETER BOX.  For a
    # point teacher this is the same number either way, but a PARAMETRIZED
    # certificate carries a y-INTERVAL, and the prefilter's dichotomy is not
    # the same law on the two sides of y = 0: positive teacher masses force
    # each student's offsets to STRADDLE the gap lattice, mixed ones force
    # them to the SAME side.  Deriving the branch from the low endpoint alone
    # would let a PREFILTER leaf be re-verified against the wrong law on a box
    # the search had rejected.  So it is derived exactly as the search derives
    # it -- from the mass interval over the whole y-box -- and a box that
    # straddles y = 0 is refused here too.
    Ybox = (J.I(str(spec["y"])) if spec.get("y_hi") is None else
            X.iv.mpf([mpmath.mpf(str(spec["y"])),
                      mpmath.mpf(str(spec["y_hi"]))]))
    s1box = X.masses_at(Ybox)[1]
    positive = float(J.lo(s1box)) > 0.0
    mixed = float(J.hi(s1box)) < 0.0
    if not (positive or mixed):
        print("SPEC REJECTED: the parameter box straddles y = 0, where the "
              "sign-law dichotomy has no single form; the search refuses such "
              "a box and so does the replay.")
        sys.exit(1)

    by_root = {}
    for key, v in leaves.items():
        ri, p = key.split(":", 1)
        by_root.setdefault(int(ri), {})[p] = v

    # The four faces of the search rectangle, regenerated from the spec by the
    # same expressions that build the root cells.  ROW 6 of the completeness
    # chain is the claim that no solution touches them.
    s_lo = mpmath.mpf(repr(spec["seam"]))
    s_hi = s_lo + J.hi(2 * X.PIv)
    d_lo = mpmath.mpf(repr(spec["delta"]))
    d_hi = J.hi(X.PIv) + mpmath.mpf(repr(spec["dmax_pad"]))

    n_ok = n_bad = n_undec = n_encl = 0
    bad = []
    margins = []
    sols = []

    def walk(ri, bx, p, depth):
        nonlocal n_ok, n_bad, n_undec, n_encl
        v = by_root.get(ri, {}).get(p)
        if v is None:
            if depth > 200:
                bad.append((ri, p, "NO LEAF, depth limit"))
                n_bad += 1
                return
            walk(ri, child(bx, 0), p + "0", depth + 1)
            walk(ri, child(bx, 1), p + "1", depth + 1)
            return
        if v.startswith("UNDECIDED"):
            n_undec += 1
            return
        ok, encl = verify_leaf(Tc, bx, v, beta_lo, beta_hi, positive)
        if ok:
            n_ok += 1
            if encl is not None:
                n_encl += 1
                S, DD = encl
                margins.append(min(J.lo(S) - s_lo, s_hi - J.hi(S),
                                   J.lo(DD) - d_lo, d_hi - J.hi(DD)))
                R1, R0 = S, S + DD
                rr = X.sep_res(Tc, R0, R1, None, Dm=DD)
                lab = (X.schur_label(Tc, R0, R1, rr[2], rr[3])[0]
                       if rr is not None else "undecided")
                sols.append({"th1_mid": float(J.mid(R1)),
                             "th0_mid": float(J.mid(R0)),
                             "D_mid": float(J.mid(DD)),
                             "width": max(J.width(R1), J.width(DD)),
                             "encl": [float(J.lo(R1)), float(J.hi(R1)),
                                      float(J.lo(DD)), float(J.hi(DD))],
                             "schur": lab})
        else:
            n_bad += 1
            if len(bad) < 20:
                bad.append((ri, p, v))

    sys.setrecursionlimit(10000)
    for ri, bx in enumerate(roots):
        walk(ri, bx, "", 0)

    print("certificate : %s" % os.path.basename(path))
    print("teacher     : beta=%s y=%s" % (spec["beta"], spec["y"]))
    print("root cells  : %d (regenerated from the spec)" % len(roots))
    print("leaves      : %d claimed, %d verified, %d FAILED, %d undecided"
          % (len(leaves), n_ok, n_bad, n_undec))
    for b in bad:
        print("   FAILED leaf root=%s path=%s verdict=%s" % b)
    if n_bad:
        print("CERTIFICATE INVALID")
        sys.exit(1)
    # ROW 6, boundary containment.  Each solution-bearing leaf carries a
    # re-derived enclosure of its solution; the margin below is the distance
    # from that enclosure to the NEAREST of the four faces of the search
    # rectangle, computed from the enclosure's outer endpoints, so it is a
    # rigorous lower bound on the distance from the true solution to the
    # boundary.  A positive minimum is exactly the hypothesis
    # `censusAngleMapJ_zeroCount_eq_witness_strip` consumes in place of the
    # schema's `hint` (four STRICT domain inequalities at every zero).
    if margins:
        mmin = min(margins)
        print("boundary    : %d solution-bearing leaves (fold and 2pi "
              "duplicates NOT merged), min distance from an enclosure to the "
              "four rectangle faces = %.6e" % (n_encl, float(mmin)))
        if mmin <= 0:
            print("BOUNDARY CONTAINMENT FAILS: a solution enclosure touches "
                  "or crosses a face of the search rectangle, so the schema's "
                  "interior hypothesis is NOT discharged here.")
            sys.exit(1)
    else:
        print("boundary    : no enclosed solution in this certificate")
    # ROW 5, the family COUNT.  The leaves above enclose solutions of the
    # (s, D) system; several leaves can enclose the SAME unordered family,
    # because the `2pi` seam and the `dmax_pad` fold each represent a pair
    # twice.  `certified_dedup` merges two enclosures only on a Krawczyk
    # `unique` verdict for their (inflated) joint hull -- a PROOF that they
    # hold one solution, never an assumed identification -- and then certifies
    # that the surviving classes are pairwise disjoint.  Re-derived here from
    # the replayed enclosures, so the count is part of the replay.
    if sols:
        d = F.certified_dedup(Tc, sols, float(spec["beta"]))
        print("families    : %d classes from %d enclosures (%d certified "
              "merges); disjointness certified: %s; exact-fit classes: %d; "
              "SEPARATED FAMILIES: %d"
              % (d["n_classes"], d["n_raw"], d["n_merges"],
                 d["disjointness_certified"], d["n_exact_fit_classes"],
                 d["n_separated_families"]))
        for c in sorted(d["classes"], key=lambda z: z["th1_mid"]):
            print("   th1=%.9f th0=%.9f D=%.9f  schur=%-10s exact_fit=%s"
                  % (c["th1_mid"], c["th0_mid"], c["D_mid"], c["schur"],
                     c["exact_fit"]))
        if not d["disjointness_certified"] or d["overlapping_class_pairs"]:
            print("FAMILY COUNT NOT CERTIFIED: classes overlap, so the count "
                  "is not a count.")
            sys.exit(1)
    if n_undec:
        print("CERTIFICATE VALID BUT INCOMPLETE: %d undecided leaves, so the "
              "domain is NOT fully resolved and no completeness claim "
              "follows." % n_undec)
        sys.exit(1)
    print("CERTIFICATE VALID AND COMPLETE: every leaf verified, none "
          "undecided, so the enclosed families are ALL of them.")


if __name__ == "__main__":
    main()
