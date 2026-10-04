"""ROUTE-B FEASIBILITY MEASUREMENT: does the PROVED sign law collapse the
per-teacher two-equation branch-and-bound that exhausted its budget in
`face_bb.py` (400000 steps at every one of the 28 face witnesses, ~95% of the
search area still undecided)?

The new pruning is PURELY ORDER-THEORETIC and costs no interval evaluation:

  sign law (Lean, SeparatedCount/SignLawJ.lean, axiom-clean)
      separatedNumJ D x < 0  for x mod 2pi in (0, D)
      separatedNumJ D x > 0  for x mod 2pi in (D, 2pi)
  student rows (SeparatedCount/ScalarColumnJ.lean)
      s0*N(D, th0)  + s1*N(D, th0-beta)  = 0
      s0*N(D, -th1) + s1*N(D, beta-th1)  = 0

  => POSITIVE teacher (s0, s1 > 0): the two offsets of each student lie on
     OPPOSITE open sides of the gap lattice {0, D} + 2pi*Z
     (`separatedPairCritical_student0/1_offsets_straddle`, LatticeStraddleJ).
  => MIXED teacher (s0 > 0 > s1): they lie on the SAME open side.

A box in (s, D) = (th1, th0-th1) whose offsets are CERTAINLY on the wrong
configuration contains no family and is discarded by four integer/endpoint
comparisons.

This script measures the collapse at ONE witness at a time and reports the
undecided-area-versus-steps curve, so the two effects are separated:
`--prefilter 0` reruns the identical search with the sign law switched OFF.
The box ordering is best-first by area (a priority queue), unlike face_bb's
LIFO stack, so the reported area curve is meaningful at any budget.

Usage:  python3 face_bb_signlaw.py <beta> <y> <budget> <0|1 prefilter> [tag]
Output: signlaw_bb_<tag>.json
"""
import json
import heapq
import math
import os
import sys
import time

import mpmath

import noncentered as J
import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
TWO_PI = 2.0 * math.pi


# --------------------------------------------------------------------------
# the sign-law prefilter (float endpoints, outward-rounded by EPS)
# --------------------------------------------------------------------------
EPS = 1e-12


def _side_certain(xa, xb, da, db):
    """Return -1 if [xa,xb] lies certainly inside (0,D) for EVERY D in
    [da,db]; +1 if certainly inside (D,2pi) for every such D; 0 otherwise.
    Uses the 2pi-shift invariance separatedNumJ_offset_int_shift."""
    if xb - xa >= da:
        return 0
    m0 = math.floor(-xa / TWO_PI)
    for m in (m0, m0 + 1, m0 + 2):
        u, v = xa + m * TWO_PI, xb + m * TWO_PI
        if u > EPS and v < da - EPS:
            return -1
        if u > db + EPS and v < TWO_PI - EPS:
            return 1
    return 0


def prefilter_kills(sa, sb, da, db, beta, positive, beta_hi=None):
    """True if the sign law excludes the whole box (s in [sa,sb],
    D in [da,db]).  `beta` may be an INTERVAL [beta, beta_hi]: the offset
    ranges are then widened outward, so a kill stays valid for EVERY teacher
    gap in that interval."""
    bl = beta
    bh = beta if beta_hi is None else beta_hi
    for (xa, xb, ya, yb) in (
            # student 0: offsets th0 = s + D and th0 - beta
            (sa + da, sb + db, sa + da - bh, sb + db - bl),
            # student 1: offsets -th1 = -s and beta - th1
            (-sb, -sa, bl - sb, bh - sa)):
        p = _side_certain(xa, xb, da, db)
        if p == 0:
            continue
        q = _side_certain(ya, yb, da, db)
        if q == 0:
            continue
        if positive and p == q:
            return True          # same side, but a straddle is forced
        if (not positive) and p != q:
            return True          # opposite sides, but same side is forced
    return False


# --------------------------------------------------------------------------
# CERTIFIED dedup of the raw Krawczyk boxes
# --------------------------------------------------------------------------
# A "unique" Krawczyk verdict on a box says: that box contains EXACTLY ONE
# solution.  Hence the merge rule, which is a proof and not a tolerance:
#
#     A unique, B unique, hull(A,B) unique   ==>   sol(A) = sol(B)
#
# (the hull contains both solutions and only one solution).  Equality is
# transitive, so union-find over this relation is sound.  Distinctness is
# the other half and needs its own certificate: two classes hold DIFFERENT
# solutions as soon as their representative enclosures are DISJOINT mod 2pi.
# When every pair is either merged or disjoint, the class count is the exact
# number of separated families of that teacher -- no visual dedup anywhere.
# The seam s = 0 ~ s = 2pi is handled by trying the +-2pi lifts of B.


def _mp(x):
    """Exact mpf from either a float or an mpf (the inflated hulls of
    `certified_same` carry mpf endpoints, the stored enclosures floats)."""
    return x if isinstance(x, mpmath.mpf) else mpmath.mpf(repr(x))


def _kraw_verdict(Tc, box):
    """Krawczyk verdict on a raw (s_lo, s_hi, d_lo, d_hi) box."""
    T1 = X.iv.mpf([_mp(box[0]), _mp(box[1])])
    Dm = X.iv.mpf([_mp(box[2]), _mp(box[3])])
    T0 = T1 + Dm
    bb = {}
    ok = True
    for key, Xi in (("n0", T0), ("n1", T1), ("n0b", T0 - Tc.B),
                    ("n1b", T1 - Tc.B), ("nD", Dm)):
        ps = J.branch_pieces(Xi)
        if len(ps) != 1:
            ok = False
            break
        bb[key] = ps[0][1]
    br = bb if ok else None
    try:
        k, _, _ = X._krawczyk_sD(Tc, T1, Dm, br)
    except Exception:                                            # noqa
        return None
    return k


def sep_res_centered_sign(Tc, T1, Dm, br):
    """Sign test on the residual pair using the MEAN-VALUE (centered) form.

    Plain interval evaluation of `(G0, G1)` loses O(w) with a large constant to
    dependency -- measured on this residual at 3.7x the true range -- so on a
    box of side w it can fail to sign a residual that is nowhere near zero.
    The mean-value form

        G(X) subset G(xc) + J(X) . (X - xc)

    is an equally valid enclosure whose dependency is second order, and it is
    available for free now that the Jacobian is trustworthy (the sign error of
    9.5 is exactly why this could not have been used before).  Returns
    `(s0, s1)` where each entry is +1 / -1 if that component is certainly of
    that sign on the box and 0 if undecided; taking the better of this and the
    plain form is sound because both are enclosures of the same quantity.

    Chain rule to (s, D) coordinates: d/ds = d/dth0 + d/dth1, d/dD = d/dth0.
    """
    T0 = T1 + Dm
    r = X.sep_res(Tc, T0, T1, br, want_jac=True, Dm=Dm)
    if r is None:
        return 0, 0
    j00, j01, j10, j11 = r[8]
    a00, a01 = j00 + j01, j00          # dG0/ds, dG0/dD
    a10, a11 = j10 + j11, j10          # dG1/ds, dG1/dD
    sm, dm = X.thin(J.mid(T1)), X.thin(J.mid(Dm))
    rc = X.sep_res(Tc, sm + dm, sm, br, Dm=dm)
    if rc is None:
        return 0, 0
    rs, rd = T1 - sm, Dm - dm
    G0 = rc[0] + a00 * rs + a01 * rd
    G1 = rc[1] + a10 * rs + a11 * rd
    return X.sgn(G0), X.sgn(G1)


def rowform_centered_sign(Tc, T1, Dm, br):
    """Mean-value sign test on the DIVISION-FREE row form.

    `E_i = massDet(D) * G_i`, so a sign-definite `E_i` excludes the box exactly
    as a sign-definite `G_i` does.  What differs is the DERIVATIVE the
    mean-value form leans on: the Cramer Jacobian divides by `massDet` three
    times per entry and then subtracts two same-sized quotients, and measured
    on 200 plateau-scale boxes at the F1 witness its `dG0/dD` is enclosed
    5041x wider than the true range (p90 2.9e5).  Differentiating `E = M*G`
    never divides; the same measurement gives 27.9x (p90 127.9), a 180x
    median improvement, and the resulting centered residual is within 5.8x of
    its true range instead of 512x.

    MEASURED PAYOFF at the F1 witness: of 400 boxes that survive the
    prefilter, the plain enclosure, the Cramer centered form AND Krawczyk,
    this test excludes 164 -- 41 % of the plateau the mixed sector was stuck
    on.  (Monotonicity on the same tightened derivative fires on only 10 of
    400: the derivative genuinely turns on these boxes, so the mean-value
    form, not monotonicity, is the right consumer of it.)

    Returns `(s0, s1)`, each +1 / -1 / 0, as `sep_res_centered_sign` does.
    """
    r = X.sep_res_rowform_jac(Tc, T1 + Dm, T1, br, Dm=Dm)
    if r is None:
        return 0, 0
    _, _, ds0, dd0, ds1, dd1 = r
    sm, dm = X.thin(J.mid(T1)), X.thin(J.mid(Dm))
    c = X.sep_res_rowform(Tc, sm + dm, sm, br, Dm=dm)
    if c is None:
        return 0, 0
    rs, rd = T1 - sm, Dm - dm
    E0 = c[0] + ds0 * rs + dd0 * rd
    E1 = c[1] + ds1 * rs + dd1 * rd
    return X.sgn(E0), X.sgn(E1)


def _certify_inflated(Tc, box, br):
    """Last-resort certification of a box that reached minimum width.

    A solution sitting exactly ON a boundary of the decomposition is never
    strictly interior to either neighbouring cell, so Krawczyk refuses both
    and they bisect down to `minw` with the family still uncertified -- the
    residual is then two adjacent boxes sharing the boundary the solution sits
    on.  Measured at the F2 witness (2.2, -0.06) and the F6 witness
    (1.5, -0.96): exactly two such boxes each, bracketing the census family at
    `D = 0.35938` resp. `0.358138`, and the certified count came out one short.

    Testing an INFLATED box is sound and it is what resolves them: if the
    inflated box `Z' ⊇ Z` contains EXACTLY ONE solution, then `Z` contains at
    most that one, so recording it loses nothing and the cell is discharged.
    Duplicates across the several boxes that share the solution are removed by
    the certified dedup, which merges only on a Krawczyk verdict.
    Returns the refined `(T1, Dm)` enclosure, or None."""
    # The ladder must start SMALL.  Krawczyk certifies only inside its
    # contraction radius, which is set by the conditioning of the Jacobian at
    # the solution: at the F2 family (D = 0.35938) the (s,D) Jacobian
    # determinant is -5.5e-3 and the test succeeds at radius 1e-5 but FAILS at
    # 1e-4 and above.  An inflation ladder starting at 1e-4 therefore never
    # fires -- measured, certified_via_inflation = 0 with the family still
    # missing.  Start below the minimum box width and climb.
    for eps in (mpmath.mpf("1e-9"), mpmath.mpf("1e-8"), mpmath.mpf("1e-7"),
                mpmath.mpf("1e-6"), mpmath.mpf("1e-5"), mpmath.mpf("1e-4")):
        z = (box[0] - eps, box[1] + eps, box[2] - eps, box[3] + eps)
        T1 = X.iv.mpf([_mp(z[0]), _mp(z[1])])
        Dm = X.iv.mpf([_mp(z[2]), _mp(z[3])])
        k, N0, N1 = X._krawczyk_sD(Tc, T1, Dm, None)
        if k != "unique":
            continue
        S, DD = N0, N1
        for _ in range(60):
            kk, M0, M1 = X._krawczyk_sD(Tc, S, DD, None)
            if kk != "unique" or (J.width(M0) >= J.width(S)
                                  and J.width(M1) >= J.width(DD)):
                break
            S, DD = M0, M1
        return S, DD
    return None


def _hull(a, b):
    return (min(a[0], b[0]), max(a[1], b[1]),
            min(a[2], b[2]), max(a[3], b[3]))


def _lift(a, k):
    return (a[0] + k * TWO_PI, a[1] + k * TWO_PI, a[2], a[3])


def _fold(a):
    """The SAME unordered student pair, written with the other student first.

    With `s = th1`, `D = th0 - th1`, the pair `{th0, th1}` is equally described
    by `s' = th0 = s + D` and `D' = 2pi - D`.  `dmax_pad` extends `D` past `pi`
    precisely so `D = pi` is interior, which makes every near-antipodal family
    appear TWICE -- once on each side of `pi` -- and the two copies are
    related by exactly this involution, not by a `2pi` shift in `s`.  Measured
    at the F5 witness `(2.6, -0.5)`: the run certified 4 disjoint boxes where
    the census has 3 families, the extra one being the fold image of the
    `D = pi` family."""
    return (a[0] + a[2], a[1] + a[3], TWO_PI - a[3], TWO_PI - a[2])


def _overlap(a, b):
    """Do the two (s, D) enclosures overlap, modulo the 2pi shift and the fold?

    Used only for the DISJOINTNESS half of the dedup, where the safe direction
    is to report an overlap: two classes count as distinct solutions only when
    no representative of one can meet the other."""
    for bb in (b, _fold(b)):
        if not (a[2] <= bb[3] and bb[2] <= a[3]):
            continue
        for k in (0, 1, -1):
            c = _lift(bb, k)
            if a[0] <= c[1] and c[0] <= a[1]:
                return True
    return False


def certified_same(Tc, a, b):
    """PROOF that two certified boxes hold the same solution.

    Tried modulo the `2pi` shift of `s` AND modulo the fold `_fold`, since the
    padded `D`-range represents each near-antipodal pair twice.  A merge is
    only ever accepted on a Krawczyk `unique` verdict for the joint hull, so
    the fold is a candidate identification, never an assumed one."""
    for bb in (b, _fold(b)):
        for k in (0, 1, -1):
            h = _hull(a, _lift(bb, k))
            if h[1] - h[0] > math.pi or h[3] - h[2] > math.pi:
                continue
            # EPSILON-INFLATION.  Krawczyk needs K(X) STRICTLY inside X, which
            # a hull of two 1e-46-wide enclosures can never satisfy -- it has
            # essentially zero width, so the strict containment fails on a box
            # that certainly holds the solution.  Testing an INFLATED hull is
            # sound and stronger: a `unique` verdict there says the inflated
            # box holds exactly one solution, and both originals lie inside
            # it, so they hold the same one.  Measured without this: the two
            # copies of the D = pi family at the F5 witness were correctly
            # flagged as overlapping but could not be merged.
            for eps in (mpmath.mpf("1e-12"), mpmath.mpf("1e-9"),
                        mpmath.mpf("1e-6"), mpmath.mpf("1e-4")):
                hi = (h[0] - eps, h[1] + eps, h[2] - eps, h[3] + eps)
                if _kraw_verdict(Tc, hi) == "unique":
                    return True
    return False


def certified_dedup(Tc, sols, beta, lat_tol=1e-5):
    """Union-find over `certified_same`, then a disjointness certificate."""
    n = len(sols)
    boxes = [tuple(z["encl"]) for z in sols]
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    merges = 0
    for i in range(n):
        for j in range(i + 1, n):
            if find(i) == find(j):
                continue
            if certified_same(Tc, boxes[i], boxes[j]):
                parent[find(i)] = find(j)
                merges += 1
    classes = {}
    for i in range(n):
        classes.setdefault(find(i), []).append(i)

    reps = []
    for root, members in classes.items():
        k = min(members, key=lambda t: (boxes[t][1] - boxes[t][0])
                + (boxes[t][3] - boxes[t][2]))
        reps.append({"members": members, "encl": list(boxes[k]),
                     "th1_mid": sols[k]["th1_mid"],
                     "th0_mid": sols[k]["th0_mid"],
                     "D_mid": sols[k]["D_mid"],
                     "width": sols[k]["width"],
                     "schur": sols[k]["schur"]})

    overlaps = []
    for i in range(len(reps)):
        for j in range(i + 1, len(reps)):
            if _overlap(tuple(reps[i]["encl"]), tuple(reps[j]["encl"])):
                overlaps.append([i, j])

    def on_lattice(lo, hi):
        """Is the whole enclosure within lat_tol of a teacher atom?"""
        for c in (0.0, beta % TWO_PI):
            for k in (-1, 0, 1):
                cc = c + k * TWO_PI
                if abs(lo - cc) < lat_tol and abs(hi - cc) < lat_tol:
                    return True
        return False

    for r in reps:
        lo1, hi1, lod, hid = r["encl"]
        r["exact_fit"] = bool(on_lattice(lo1, hi1)
                              and on_lattice(lo1 + lod, hi1 + hid))
    n_fit = sum(1 for r in reps if r["exact_fit"])
    return {
        "n_raw": n, "n_merges": merges, "n_classes": len(reps),
        "classes": reps,
        "overlapping_class_pairs": overlaps,
        "disjointness_certified": not overlaps,
        "n_exact_fit_classes": n_fit,
        "n_separated_families": len(reps) - n_fit,
    }


# --------------------------------------------------------------------------
# the branch and bound, best-first by area
# --------------------------------------------------------------------------
def run(beta_f, y_f, budget, use_prefilter, delta=mpmath.mpf("0.02"),
        minw=mpmath.mpf("1e-8"), cell=0.35, beta_hi=None, y_hi=None,
        seam=mpmath.mpf("0.137"), dmax_pad=mpmath.mpf("0.05"),
        atom_eps=mpmath.mpf("1e-3"), order="area", emit_cert=False):
    """Enumerate the separated families of a teacher, or of a whole BOX of
    teachers when `beta_hi` / `y_hi` are given.

    PARAMETRIZED MODE is what turns a witness result into a FACE result.
    Every interval evaluation here is monotone in the parameter box, so:
      * an exclusion (`sgn(G) != 0` on the box) holds for EVERY teacher in
        the box, and
      * a Krawczyk verdict `K(Z) subset int(Z)` on a box `Z` of student
        coordinates, evaluated with the teacher entering as an interval,
        says that for EVERY teacher in the parameter box the system has
        exactly one solution in `Z`.
    So when the search finishes with an empty undecided list, the family
    COUNT is certified constant across the whole parameter box -- which is
    the count-constancy obligation, discharged by interval arithmetic rather
    than by an implicit-function theorem.  (The solutions move with the
    teacher; the enclosure covers all of them.)"""
    # DOMAIN (2026-08-08).  The obvious domain `s in [0, 2pi]`, `D in
    # [delta, pi]` puts genuine solutions ON its boundary: the EXACT FIT sits
    # at the seam `s = 0`, and at several witnesses a family sits at exactly
    # `D = pi` (the antipodal stratum).  A solution on the boundary can never
    # be Krawczyk-certified -- the test needs it strictly interior -- and the
    # boxes touching it can never be excluded either, so they subdivide
    # forever and the search cannot terminate.  Two changes fix that:
    #   `seam`      shifts the s-period to [seam, seam + 2pi] so the exact fit
    #               is interior (any generic offset does);
    #   `dmax_pad`  extends D past pi to pi + pad, using the fold
    #               D <-> 2pi - D.  Pairs with D in (pi, pi + pad) are mirror
    #               duplicates of pairs just below pi, so this OVER-counts
    #               rather than misses -- the safe direction for a
    #               completeness claim -- and it makes D = pi interior.
    B = J.I(str(beta_f)) if beta_hi is None else X.iv.mpf(
        [mpmath.mpf(str(beta_f)), mpmath.mpf(str(beta_hi))])
    Y = J.I(str(y_f)) if y_hi is None else X.iv.mpf(
        [mpmath.mpf(str(y_f)), mpmath.mpf(str(y_hi))])
    s0, s1 = X.masses_at(Y)
    Tc = X.Teacher(B, s0, s1)
    positive = float(J.lo(s1)) > 0.0
    mixed = float(J.hi(s1)) < 0.0
    if not (positive or mixed):
        raise SystemExit("parameter box straddles y = 0: the sign-law "
                         "dichotomy has no single form on it; split the box")
    beta = float(J.lo(B))
    beta_h = float(J.hi(B))

    two_pi = J.hi(2 * X.PIv)
    dpi = J.hi(X.PIv)
    s_lo, s_hi = seam, seam + two_pi
    d_hi = dpi + dmax_pad
    cuts = {s_lo, s_hi}
    a, b = J.endpoints(Tc.B)
    for z in (a, b, a + dpi, b + dpi, a - dpi, b - dpi,
              a + two_pi, b + two_pi, dpi + seam):
        if s_lo < z < s_hi:
            cuts.add(z)
    # The cut set collects teacher atom lines and their half/full-turn
    # translates, several of which can coincide to within rounding.  A
    # duplicated cut makes a cell of ZERO width, which has zero area, can
    # never be excluded, and (under depth-first) is bisected forever in the
    # other coordinate -- measured: 691 minwidth records at area 5.9e-54, all
    # on the single line s = beta.  Dedupe with a tolerance and keep only
    # strictly positive cells; the atom lines then sit on cell BOUNDARIES,
    # which is what the branch bookkeeping wants anyway.
    def _dedupe(cs, tol=mpmath.mpf("1e-12")):
        cs = sorted(cs)
        out = [cs[0]]
        for z in cs[1:]:
            if z - out[-1] > tol:
                out.append(z)
        return out

    # ATOM LINES MUST NOT BE CELL BOUNDARIES.  Same defect as pi was: a
    # solution lying exactly ON a cut can never be Krawczyk-certified, and the
    # cells on either side of it can never be excluded.  Measured at the F5
    # witness: after every other fix the last two undecided boxes sat at
    # s = 0 mod 2pi, D = beta -- the EXACT FIT, both students on the teacher
    # lattice -- straddling the cut D = beta.  So each atom-derived cut `z` is
    # replaced by the PAIR `z - atom_eps, z + atom_eps`, putting the atom line
    # strictly inside a thin cell.  The thin cell straddles the kernel's kink,
    # where `branch_pieces` gives more than one piece and `sep_res` falls back
    # to the branch-free atoms -- a wider but still valid enclosure.
    cuts = {z for z in cuts if z in (s_lo, s_hi)} | {
        w for z in cuts if z not in (s_lo, s_hi)
        for w in (z - atom_eps, z + atom_eps) if s_lo < w < s_hi}
    cuts.add(s_lo)
    cuts.add(s_hi)
    cuts = _dedupe(cuts)
    sf = []
    for i in range(len(cuts) - 1):
        u, v = cuts[i], cuts[i + 1]
        n = max(1, int(math.ceil(float(v - u) / cell)))
        for k in range(n):
            lo, hi = u + (v - u) * k / n, u + (v - u) * (k + 1) / n
            if hi > lo:
                sf.append((lo, hi))
    # NOTE: `dpi` is deliberately NOT a cut.  Putting pi on a cell BOUNDARY
    # is exactly what dmax_pad was introduced to avoid -- the antipodal
    # family sits at D = pi exactly, and a solution on a box boundary can
    # never be Krawczyk-certified.  Measured before this line was fixed:
    # depth-first drilled both F4 solutions to the 2e-5 minimum width with D
    # pinned to [pi, pi + 1.2e-5] and no verdict.
    dcuts = list({delta, d_hi} | {z for z in (a, b, a + dpi, b + dpi,
                                              two_pi - b, two_pi - a)
                                  if delta < z < d_hi})
    dcuts = [z for z in dcuts if z in (delta, d_hi)] + [
        w for z in dcuts if z not in (delta, d_hi)
        for w in (z - atom_eps, z + atom_eps) if delta < w < d_hi]
    dcuts = _dedupe(dcuts)
    df = []
    for i in range(len(dcuts) - 1):
        u, v = dcuts[i], dcuts[i + 1]
        n = max(1, int(math.ceil(float(v - u) / cell)))
        for k in range(n):
            lo, hi = u + (v - u) * k / n, u + (v - u) * (k + 1) / n
            if hi > lo:
                df.append((lo, hi))

    def area(bx):
        return float((bx[1] - bx[0]) * (bx[3] - bx[2]))

    # ORDERING (measured 2026-08-08).  Best-first by AREA is the wrong
    # priority for TERMINATION.  A box near a true solution can never be
    # excluded (the residual really does vanish there) and can only be
    # resolved by Krawczyk, which needs a box of about the certification
    # radius -- here ~1e-3.  Area-first refines the WHOLE domain to that
    # scale before it ever reaches one, costing ~20/(1e-3)^2 = 2e7 boxes;
    # depth-first drills each surviving box straight to resolution and pays
    # ~2^18 only in the neighbourhoods that need it.  Area-first is still the
    # right choice for MEASURING the collapse (its area curve is meaningful
    # at any budget), so both are kept.
    # HYBRID (measured 2026-08-08).  The two pure orderings fail in opposite
    # ways.  Area-first covers the domain but must refine EVERYWHERE down to
    # the certification scale before it can certify anything, which is
    # hopeless once minw is small.  Depth-first resolves a neighbourhood
    # immediately but drills every surviving box to full depth before looking
    # at the next one, so on the mixed-sign faces -- where the order-only
    # prefilter removes only ~39% and the rest must be excluded by the
    # residual -- coverage crawls: measured at the F1 witness, live area
    # 12.79 of 19.9 still standing after 105000 depth-first steps.
    #
    # The hybrid takes SMALL boxes first (they are cheap and each one either
    # excludes, certifies, or drills a solution neighbourhood to resolution)
    # and otherwise the LARGEST box (coverage).  So resolution happens at
    # depth-first speed while the frontier still advances area-first.
    depth_first = (order == "depth")
    hybrid = (order == "hybrid")
    dive_w = mpmath.mpf("1e-3")

    def prio(bx):
        if hybrid:
            small = (bx[1] - bx[0]) < dive_w and (bx[3] - bx[2]) < dive_w
            return (0, float(area(bx))) if small else (1, -float(area(bx)))
        if depth_first:
            return (0, float(area(bx)))
        return (0, -float(area(bx)))
    # THE CERTIFICATE (owner artifact standard, 2026-08-08).  The SEARCH --
    # the heap, the priorities, the orderings -- is prospecting and no number
    # from it may stand in the completeness chain.  What is emitted instead is
    # a replayable object: the bisection TREE, as the set of leaf PATHS from
    # each root cell together with each leaf's verdict.  Paths rather than
    # coordinates, because the split rule is deterministic (bisect the wider
    # side), so the checker regenerates every box exactly and needs no float
    # comparison to identify it.  Verifying it needs none of the search: walk
    # the tree, and at each leaf re-derive the claimed verdict.  Coverage of
    # the domain is then structural -- a node is a leaf or the union of its
    # two children -- so a tree whose every leaf carries a valid verdict IS
    # the completeness proof.
    leaves = {} if emit_cert else None
    heap = []
    seq = 0
    total_area = 0.0
    roots = []
    for x in sf:
        for z in df:
            bx = (x[0], x[1], z[0], z[1])
            roots.append(bx)
            total_area += area(bx)
            p0, p1 = prio(bx)
            heapq.heappush(heap, (p0, p1, seq, (bx, len(roots) - 1, "")))
            seq += 1

    sols = []
    undecided = []
    steps = 0
    n_excl_plain = 0
    n_inflated = 0
    n_excl_centered = 0
    killed_prefilter = 0
    area_prefilter = 0.0
    curve = []
    t_start = time.time()

    def live_area():
        return (sum(area(k[3][0]) for k in heap)
                + sum(u["area"] for u in undecided))

    while heap:
        steps += 1
        if steps > budget:
            for k in heap:
                if leaves is not None:
                    leaves["%d:%s" % (k[3][1], k[3][2])] = "UNDECIDED_BUDGET"
                undecided.append({"box": [float(t) for t in k[3][0]],
                                  "area": area(k[3][0]),
                                  "reason": "budget"})
            heap = []
            break
        if steps % 5000 == 1:
            curve.append({"steps": steps, "live_area": live_area(),
                          "heap": len(heap), "undecided": len(undecided),
                          "sols": len(sols),
                          "seconds": time.time() - t_start})
            print("  step %7d  live_area %9.5f  heap %6d  undec %5d  "
                  "sols %4d  %6.0fs"
                  % (steps, curve[-1]["live_area"], len(heap),
                     len(undecided), len(sols), curve[-1]["seconds"]),
                  flush=True)
        _, _, _, (bx, root_i, path) = heapq.heappop(heap)
        s0_, s1_, d0_, d1_ = bx

        def _leaf(tag):
            if leaves is not None:
                leaves["%d:%s" % (root_i, path)] = tag

        if use_prefilter and prefilter_kills(float(s0_), float(s1_),
                                             float(d0_), float(d1_),
                                             beta, positive, beta_h):
            killed_prefilter += 1
            area_prefilter += area(bx)
            _leaf("PREFILTER")
            continue

        T1 = X.iv.mpf([s0_, s1_])
        Dm = X.iv.mpf([d0_, d1_])
        T0 = T1 + Dm
        bb = {}
        ok = True
        for key, Xi in (("n0", T0), ("n1", T1), ("n0b", T0 - Tc.B),
                        ("n1b", T1 - Tc.B), ("nD", Dm)):
            ps = J.branch_pieces(Xi)
            if len(ps) != 1:
                ok = False
                break
            bb[key] = ps[0][1]
        br = bb if ok else None
        r = X.sep_res(Tc, T0, T1, br, Dm=Dm)
        split = False
        if r is None:
            if (s1_ - s0_) < minw and (d1_ - d0_) < minw:
                _leaf("UNDECIDED_MASSDET")
                undecided.append({"box": [float(t) for t in bx],
                                  "area": area(bx),
                                  "reason": "massDet_zero"})
                continue
            split = True
        elif X.sgn(r[0]) != 0 or X.sgn(r[1]) != 0:
            n_excl_plain += 1
            _leaf("EXCL_PLAIN")
            continue
        # ORDER: the ROW-FORM centered test first.  Measured on 400 boxes of
        # this face, it strictly DOMINATES the Cramer one -- 317 boxes killed
        # by both, 21 by the row form alone, ZERO by the Cramer form alone --
        # and it is also marginally cheaper (4.4 ms against 4.6 ms).  The
        # Cramer test is kept as a fallback rather than deleted because
        # dominance is measured, not proved; in that role it costs only the
        # ~15 % of boxes where the row form fails.
        elif (any(rowform_centered_sign(Tc, T1, Dm, br))
              or any(sep_res_centered_sign(Tc, T1, Dm, br))):
            # the plain enclosure could not sign either component, the
            # centered one can: no solution in this box
            n_excl_centered += 1
            _leaf("EXCL_CENTERED")
            continue
        else:
            k, N0, N1 = X._krawczyk_sD(Tc, T1, Dm, br)
            if k == "empty":
                _leaf("KRAWCZYK_EMPTY")
                continue
            if k == "unique":
                S, DD = N0, N1
                for _ in range(60):
                    kk, M0, M1 = X._krawczyk_sD(Tc, S, DD, br)
                    if kk != "unique" or (J.width(M0) >= J.width(S)
                                          and J.width(M1) >= J.width(DD)):
                        break
                    S, DD = M0, M1
                R1 = S
                R0 = S + DD
                rr = X.sep_res(Tc, R0, R1, br, Dm=DD)
                lab, sg = X.schur_label(Tc, R0, R1, rr[2], rr[3])
                _leaf("KRAWCZYK_UNIQUE")
                sols.append({"th1_mid": float(J.mid(R1)),
                             "th0_mid": float(J.mid(R0)),
                             "D_mid": float(J.mid(DD)),
                             "width": max(J.width(R1), J.width(DD)),
                             "encl": [float(J.lo(R1)), float(J.hi(R1)),
                                      float(J.lo(DD)), float(J.hi(DD))],
                             "seed_box": [float(t) for t in bx],
                             "path": "%d:%s" % (root_i, path),
                             "schur": lab,
                             "sgn_c": [X.sgn(rr[2]), X.sgn(rr[3])]})
                continue
            if (s1_ - s0_) < minw and (d1_ - d0_) < minw:
                inf = _certify_inflated(Tc, bx, br)
                if inf is not None:
                    S, DD = inf
                    R1, R0 = S, S + DD
                    rr = X.sep_res(Tc, R0, R1, None, Dm=DD)
                    lab, sg = X.schur_label(Tc, R0, R1, rr[2], rr[3])
                    n_inflated += 1
                    _leaf("KRAWCZYK_INFLATED")
                    sols.append({"th1_mid": float(J.mid(R1)),
                                 "th0_mid": float(J.mid(R0)),
                                 "D_mid": float(J.mid(DD)),
                                 "width": max(J.width(R1), J.width(DD)),
                                 "encl": [float(J.lo(R1)), float(J.hi(R1)),
                                          float(J.lo(DD)), float(J.hi(DD))],
                                 "seed_box": [float(t) for t in bx],
                                 "schur": lab,
                                 "via": "inflated",
                                 "path": "%d:%s" % (root_i, path),
                                 "sgn_c": [X.sgn(rr[2]), X.sgn(rr[3])]})
                    continue
                _leaf("UNDECIDED_MINWIDTH")
                undecided.append({"box": [float(t) for t in bx],
                                  "area": area(bx), "reason": "minwidth"})
                continue
            split = True
        if split:
            if s1_ - s0_ >= d1_ - d0_:
                m = (s0_ + s1_) / 2
                kids = ((s0_, m, d0_, d1_), (m, s1_, d0_, d1_))
            else:
                m = (d0_ + d1_) / 2
                kids = ((s0_, s1_, d0_, m), (s0_, s1_, m, d1_))
            for bit, nb in enumerate(kids):
                p0, p1 = prio(nb)
                heapq.heappush(heap, (p0, p1, seq,
                                      (nb, root_i, path + str(bit))))
                seq += 1

    dedup = certified_dedup(Tc, sols, beta) if sols else None

    return {
        "beta": beta_f, "y": y_f, "budget": budget,
        "beta_hi": beta_hi, "y_hi": y_hi,
        "order": order,
        "n_root_cells": len(roots),
        "cert_leaves": leaves,
        "cell": float(cell),
        "minw": float(minw),
        "delta": float(delta), "seam": float(seam),
        "atom_eps": float(atom_eps),
        "dmax_pad": float(dmax_pad),
        "parametrized": bool(beta_hi is not None or y_hi is not None),
        "dedup": dedup,
        "prefilter": bool(use_prefilter),
        "positive_teacher": bool(positive),
        "total_area": total_area,
        "steps": steps,
        "complete": bool(steps <= budget),
        "n_sols": len(sols), "sols": sols,
        "n_undecided": len(undecided),
        "undecided": undecided,
        "undecided_area": sum(u["area"] for u in undecided),
        "undecided_reasons": {rr: sum(1 for u in undecided
                                      if u["reason"] == rr)
                              for rr in ("budget", "minwidth",
                                         "massDet_zero")},
        "certified_via_inflation": n_inflated,
        "excluded_plain_form": n_excl_plain,
        "excluded_centered_form": n_excl_centered,
        "prefilter_kills": killed_prefilter,
        "prefilter_area": area_prefilter,
        "curve": curve,
        "seconds": time.time() - t_start,
    }


if __name__ == "__main__":
    bf = float(sys.argv[1])
    yf = float(sys.argv[2])
    bud = int(sys.argv[3])
    pf = int(sys.argv[4])
    tag = sys.argv[5] if len(sys.argv) > 5 else ("%s_%s_%d" % (bf, yf, pf))
    # optional parameter BOX: beta in [bf, bhi], y in [yf, yhi]
    bhi = float(sys.argv[6]) if len(sys.argv) > 6 and sys.argv[6] != "-" \
        else None
    yhi = float(sys.argv[7]) if len(sys.argv) > 7 and sys.argv[7] != "-" \
        else None
    # the coincidence collar: families with D below it are OUT of the domain,
    # so it must sit below the smallest gap the face carries (at beta = 3.10
    # every family has D <= 0.064, so 0.02 was excluding real ones)
    dlt = mpmath.mpf(sys.argv[8]) if len(sys.argv) > 8 \
        else mpmath.mpf("0.02")
    print("beta=%s%s y=%s%s budget=%d prefilter=%d"
          % (bf, "" if bhi is None else "..%s" % bhi,
             yf, "" if yhi is None else "..%s" % yhi, bud, pf), flush=True)
    ordr = sys.argv[9] if len(sys.argv) > 9 else "area"
    out = run(bf, yf, bud, pf, delta=dlt, beta_hi=bhi, y_hi=yhi, order=ordr,
              emit_cert=True)
    cert = out.pop("cert_leaves", None)
    with open(os.path.join(HERE, "signlaw_bb_%s.json" % tag), "w") as fh:
        json.dump(out, fh, indent=1)
    if cert is not None:
        spec = {k: out[k] for k in
                ("beta", "y", "beta_hi", "y_hi", "delta", "seam", "dmax_pad",
                 "atom_eps", "minw", "cell", "n_root_cells")}
        spec["n_leaves"] = len(cert)
        with open(os.path.join(HERE, "ivcert_%s.json" % tag), "w") as fh:
            json.dump({"spec": spec, "leaves": cert}, fh, separators=(",", ":"))
        print("wrote ivcert_%s.json  (%d leaves)" % (tag, len(cert)))
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("curve", "sols", "dedup",
                                   "undecided")}, indent=1))
    d = out["dedup"]
    if d:
        print("CERTIFIED DEDUP: %d raw boxes -> %d classes (%d merges); "
              "disjointness certified: %s; exact-fit classes: %d; "
              "SEPARATED FAMILIES: %d"
              % (d["n_raw"], d["n_classes"], d["n_merges"],
                 d["disjointness_certified"], d["n_exact_fit_classes"],
                 d["n_separated_families"]))
        for c in sorted(d["classes"], key=lambda z: z["th1_mid"]):
            print("   th1=%.9f th0=%.9f D=%.9f  w=%.2e  schur=%-18s "
                  "exact_fit=%s  (%d raw boxes)"
                  % (c["th1_mid"], c["th0_mid"], c["D_mid"], c["width"],
                     c["schur"], c["exact_fit"], len(c["members"])))
        if d["overlapping_class_pairs"]:
            print("   OVERLAPPING (not separated, not merged):",
                  d["overlapping_class_pairs"])
    print("UNDECIDED BOXES:", out["n_undecided"], "area",
          out["undecided_area"])
    for u in out["undecided"][:12]:
        print("   ", u["reason"], [round(t, 9) for t in u["box"]],
              "area %.2e" % u["area"])
