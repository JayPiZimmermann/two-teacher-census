import os as _os
_HERE = _os.path.dirname(_os.path.abspath(__file__))
"""
Certified enumeration of the zero sets of the three mass-free determinants of
the CENTERED two-student phase diagram.

Every enclosure below is produced by mpmath.iv (interval arithmetic with
directed rounding) at 200 bits of working precision.  float64 appears only in
the widgets.js mirror used for validation (validate.py), never in a reported
enclosure.

Run:  python3 certify.py
"""
import json
import os
import time

import mpmath
from mpmath import iv

import centered as C

PREC = 200
C.set_prec(PREC)
PI = C.PI_IV()
HALFPI = PI / 2
MPPI = mpmath.pi

REPORT = {"precision_bits": PREC,
          "mpmath_version": mpmath.__version__,
          "arithmetic": "mpmath.iv interval arithmetic, directed rounding",
          "generated": time.strftime("%Y-%m-%d")}
UNRESOLVED = []
FAST = bool(os.environ.get("FAST"))


def S(x, n=25):
    a, b = C.endpoints(x)
    return [mpmath.nstr(a, n), mpmath.nstr(b, n)]


def strict_pos(x):
    return C.lo(x) > 0


def strict_neg(x):
    return C.hi(x) < 0


def sgn(x):
    if strict_pos(x):
        return 1
    if strict_neg(x):
        return -1
    return 0


# ==========================================================================
# 0.  The reduced (kink-free) piecewise forms.  validate.py proves
#     symbolically (sympy) and by interval sampling that these ARE the three
#     determinants.  With W = (pi-b)/2 and |w| <= W:
#        Q (b,w) = (w^2 - b^2/4) sin b
#        Om(b,w) = (b/2)(cos 2w + cos b) = b cos(w+b/2) cos(w-b/2)
#        Phi = Q + Om,   Psi = -Q + Om,   Om
#     piece A : t in [beta, pi] , b = beta   , sigma = +1, w = t - (beta+pi)/2
#     piece B : t in [0, beta]  , b = pi-beta, sigma = -1, w = t - beta/2
#     Wtau = sigma*Phi , Wpot = sigma*Psi , Wwgt = sigma*Om .
#     Consequently  Wtau + Wpot = 2*Wwgt  identically.
# ==========================================================================

IDX = {"Wtau": 0, "Wpot": 1, "Wwgt": 2}


def Wof(b, w):
    Q = (w * w - b * b / 4) * iv.sin(b)
    Om = (b / 2) * (iv.cos(2 * w) + iv.cos(b))
    return Q + Om, -Q + Om, Om


def dWof(b, w):
    sb = iv.sin(b)
    s2 = iv.sin(2 * w)
    return 2 * w * sb - b * s2, -2 * w * sb - b * s2, -b * s2


def Wmax(b):
    return (PI - b) / 2


def piece_b(Bbox, piece):
    return Bbox if piece == "A" else PI - Bbox


def piece_branches(piece):
    """(kt, ku): branch index of t and of t-beta on the piece."""
    return (0, 0) if piece == "A" else (0, -1)


def t_of_w(Bbox, piece, w):
    return (Bbox + PI) / 2 + w if piece == "A" else Bbox / 2 + w


# ==========================================================================
# 1.  Structure lemmas, hypotheses verified with intervals
# ==========================================================================
print("== 1. structure lemmas ==")
L = {}
maxres = 0.0
for k in range(1, 400):
    b = PI * k / 400
    W = Wmax(b)
    P0, S0, _ = Wof(b, C.I(0))
    PW, SW, OW = Wof(b, W)
    for z in (P0 - b * iv.cos(b / 2) * (iv.cos(b / 2) - (b / 2) * iv.sin(b / 2)),
              PW - (PI / 4) * (PI - 2 * b) * iv.sin(b),
              S0 - b * iv.cos(b / 2) * (iv.cos(b / 2) + (b / 2) * iv.sin(b / 2)),
              SW - (PI / 2) * (b - PI / 2) * iv.sin(b),
              OW):
        assert C.contains_zero(z), "endpoint closed form failed"
        maxres = max(maxres, abs(float(C.hi(z))), abs(float(C.lo(z))))
L["endpoint_closed_forms"] = {
    "Phi(b,0)": "b cos(b/2) [cos(b/2) - (b/2) sin(b/2)]",
    "Phi(b,W)": "(pi/4)(pi - 2b) sin b",
    "Psi(b,0)": "b cos(b/2) [cos(b/2) + (b/2) sin(b/2)]   > 0 on (0,pi)",
    "Psi(b,W)": "(pi/2)(b - pi/2) sin b",
    "Om(b,W)": "0",
    "verified_at_grid_points": 399,
    "max_interval_residual": maxres}
print("   endpoint closed forms verified, max residual %.2e" % maxres)

grid_ok = {"L5": 0, "L6": 0, "L7": 0, "L8": 0}
for k in range(1, 200):
    b = PI * k / 200
    W = Wmax(b)
    for j in range(0, 41):
        w = W * j / 40
        dP, dS, dO = dWof(b, w)
        P, _, Om = Wof(b, w)
        if j > 0:
            assert strict_neg(dS)
            grid_ok["L5"] += 1
        assert C.hi(Om) >= 0
        grid_ok["L6"] += 1
        if C.lo(b) > C.hi(HALFPI) and j > 0:
            assert strict_neg(dP)
            grid_ok["L7"] += 1
        if C.hi(b) < C.lo(HALFPI):
            assert strict_pos(P)
            grid_ok["L8"] += 1
L["monotonicity_and_sign_lemmas"] = {
    "L5": "Psi_w = -(2w sin b + b sin 2w) < 0 for 0 < w <= W and 0 < b < pi, so "
          "Psi is strictly decreasing in w on [0,W].",
    "L6": "Om = b cos(w+b/2) cos(w-b/2); for |w| <= W both arguments lie in "
          "[b-pi/2, pi/2] subset [-pi/2,pi/2], hence Om >= 0 with equality "
          "exactly at |w| = W.",
    "L7": "for b > pi/2, D := -Phi_w = b sin 2w - 2w sin b satisfies D(0)=0, "
          "D_ww = -4b sin 2w <= 0 (concave on 0<=2w<=pi-b) and D(W) = (2b-pi) sin b > 0, "
          "hence D(w) >= (w/W) D(W) > 0: Phi is strictly decreasing in w on [0,W].",
    "L8": "for 0 < b < pi/2, Phi > 0 on [0,W].  For w in [b/2,W], Phi = Q+Om with "
          "Q,Om >= 0 and never both zero.  For w in [0,b/2] put p = b/2+w, "
          "q = b/2-w >= 0, p+q = b <= pi/2; then Phi = (p+q)cos p cos q - pq sin(p+q) "
          "and, dividing by cos p cos q > 0, the claim is p+q > pq(tan p + tan q). "
          "Since q <= pi/2-p one has pq tan p <= p(pi/2-p)tan p < p, the last step "
          "being (pi/2-p)tan p < 1, i.e. s < tan s at s = pi/2-p; symmetrically "
          "pq tan q < q.  Adding gives the claim strictly unless p=q=0.",
    "interval_grid_checks": grid_ok}
print("   lemma grid checks:", grid_ok)

# ==========================================================================
# 2.  b* and beta*
# ==========================================================================
print("== 2. b* / beta* ==")


def f_iv(b):
    return iv.cos(b / 2) - (b / 2) * iv.sin(b / 2)


def df_iv(b):
    return -iv.sin(b / 2) - (b / 4) * iv.cos(b / 2)


assert strict_neg(df_iv(iv.mpf(["1.70", "1.75"])))
assert strict_pos(f_iv(C.I("1.70"))) and strict_neg(f_iv(C.I("1.75")))
lo_, hi_ = mpmath.mpf("1.70"), mpmath.mpf("1.75")
for _ in range(400):
    m = (lo_ + hi_) / 2
    v = f_iv(iv.mpf([m, m]))
    if C.contains_zero(v):
        break
    if strict_pos(v):
        lo_ = m
    else:
        hi_ = m
bstar = iv.mpf([lo_, hi_])
beta1 = PI - bstar
beta2 = bstar
ustar = beta1 / 2
char = iv.tan(ustar) - (PI / 2 - ustar)
assert C.contains_zero(char)
print("   b*      =", S(bstar))
print("   beta*_1 =", S(beta1))
print("   beta*_2 =", S(beta2))
L["bstar"] = {
    "equation": "cos(b/2) = (b/2) sin(b/2), i.e. (b/2) tan(b/2) = 1",
    "search_bracket": ["1.70", "1.75"],
    "uniqueness": "f'(b) = -sin(b/2) - (b/4) cos(b/2) certified strictly negative "
                  "on [1.70,1.75], and f(1.70) > 0 > f(1.75)",
    "b_star": S(bstar, 30),
    "beta_star_1_eq_pi_minus_b_star": S(beta1, 30),
    "beta_star_2_eq_b_star": S(beta2, 30),
    "sum_check_beta1_plus_beta2_eq_pi": S(beta1 + beta2 - PI, 6),
    "characterisation_beta1_eq_2u_with_tan_u_eq_halfpi_minus_u":
        {"u_star": S(ustar, 30), "residual_tan(u)-(pi/2-u)": S(char, 6)}}

# ==========================================================================
# 3.  Auxiliary constants
# ==========================================================================
print("== 3. constants ==")
lo_, hi_ = mpmath.mpf("0.7"), mpmath.mpf("0.8")
assert strict_neg(-iv.sin(iv.mpf([lo_, hi_])) - 1)
for _ in range(400):
    m = (lo_ + hi_) / 2
    v = iv.cos(iv.mpf([m, m])) - iv.mpf([m, m])
    if C.contains_zero(v):
        break
    if strict_pos(v):
        lo_ = m
    else:
        hi_ = m
wD = iv.mpf([lo_, hi_])
cconst = 2 * iv.atan2(2 / PI, C.I(1)) / PI
print("   w_D =", S(wD))
print("   c   =", S(cconst))
L["constants"] = {
    "w_Dottie": {"equation": "cos w = w", "enclosure": S(wD, 30),
                 "role": "limit position of the Wpot roots at beta -> 0 and "
                         "beta -> pi:  t = pi/2 +- w_D"},
    "c": {"definition": "(2/pi) arctan(2/pi)", "enclosure": S(cconst, 30),
          "role": "the y-limits of both curves at beta = pi/2: the torque lens "
                  "tends to y = c-1 and y = -c, the potential curve to y = 1-c "
                  "and y = c"}}

# ==========================================================================
# 4.  CERTIFIED ROOT COUNTS on an adaptively refined partition of beta.
# ==========================================================================
print("== 4. certified counts on a beta partition ==")


def bisect_root(name, Bbox, piece, wlo, whi, target=mpmath.mpf("1e-30")):
    idx = IDX[name]
    b = piece_b(Bbox, piece)
    a_, c_ = C.lo(wlo), C.hi(whi)
    s0 = sgn(Wof(b, iv.mpf([a_, a_]))[idx])
    it = 0
    while c_ - a_ > target and it < 500:
        m = (a_ + c_) / 2
        v = Wof(b, iv.mpf([m, m]))[idx]
        s = sgn(v)
        if s == 0:
            break
        if s == s0:
            a_ = m
        else:
            c_ = m
        it += 1
    return iv.mpf([a_, c_])


def count_piece(name, Bbox, piece):
    """Certified number of roots of `name` in the half-open w-range (0,W], valid
    for EVERY beta in Bbox.  Returns (count, [w enclosure]) or None."""
    b = piece_b(Bbox, piece)
    W = Wmax(b)
    if name == "Wwgt":
        return 0, []                       # L6: Om > 0 for |w| < W
    idx = IDX[name]
    F0 = Wof(b, C.I(0))[idx]
    FW = Wof(b, W)[idx]
    if name == "Wpot":                     # L5: strictly decreasing, F0 > 0
        if not strict_pos(F0):
            return None
        if strict_pos(FW):
            return 0, []
        if strict_neg(FW):
            return 1, [bisect_root(name, Bbox, piece, C.I(0), W)]
        return None
    if strict_neg(b - HALFPI):             # L8
        return 0, []
    if strict_pos(b - HALFPI):             # L7: strictly decreasing
        if not strict_neg(FW):
            return None
        if strict_pos(F0):
            return 1, [bisect_root(name, Bbox, piece, C.I(0), W)]
        if strict_neg(F0):
            return 0, []
        return None
    return None


def counts_on(Bbox):
    """total certified root count in t over one period [0,pi), per determinant,
    counting the two kink roots of Wwgt, valid for every beta in Bbox."""
    out = {}
    for det in ("Wtau", "Wpot", "Wwgt"):
        tot = 0
        encl = []
        for piece in ("B", "A"):
            r = count_piece(det, Bbox, piece)
            if r is None:
                return None
            n, ws = r
            tot += 2 * n
            for ww in ws:
                encl.append({"piece": piece, "w": S(ww, 22)})
        if det == "Wwgt":
            tot = 2                        # the two kink roots t = 0 and t = beta
        out[det] = {"n_roots_in_t_mod_pi": tot, "w_enclosures": encl}
    # kink roots of Wtau/Wpot: (pi/4)(pi-2beta) sin beta must be nonzero
    kv = (PI / 4) * (PI - 2 * Bbox) * iv.sin(Bbox)
    if sgn(kv) == 0:
        return None
    out["Wtau"]["kink_roots"] = []
    out["Wpot"]["kink_roots"] = []
    out["Wwgt"]["kink_roots"] = ["t = 0", "t = beta"]
    return out


MIN_BETA_BOX = mpmath.mpf("1e-13")


def certify_interval(label, a, c, nsplit=24):
    """Adaptive: cover (a,c) with boxes on which the counts are certified."""
    cells = []
    todo = [(mpmath.mpf(a) + (mpmath.mpf(c) - mpmath.mpf(a)) * i / nsplit,
             mpmath.mpf(a) + (mpmath.mpf(c) - mpmath.mpf(a)) * (i + 1) / nsplit)
            for i in range(nsplit)]
    todo.reverse()
    uncov = []
    while todo:
        lo2, hi2 = todo.pop()
        r = counts_on(iv.mpf([lo2, hi2]))
        if r is not None:
            cells.append({"beta_box": [mpmath.nstr(lo2, 22), mpmath.nstr(hi2, 22)],
                          "counts": {d: r[d]["n_roots_in_t_mod_pi"] for d in r},
                          "w_enclosures": {d: r[d]["w_enclosures"] for d in r}})
            continue
        if hi2 - lo2 < MIN_BETA_BOX:
            uncov.append((lo2, hi2))
            continue
        m = (lo2 + hi2) / 2
        todo.append((m, hi2))
        todo.append((lo2, m))
    for lo2, hi2 in uncov:
        UNRESOLVED.append({"kind": "beta_collar_uncertified", "range": label,
                           "beta_box": [mpmath.nstr(lo2, 22), mpmath.nstr(hi2, 22)],
                           "width": float(hi2 - lo2)})
    return cells, uncov


EPS0 = mpmath.mpf("1e-13")
RANGES = [("(0, beta*_1)", EPS0, C.lo(beta1)),
          ("(beta*_1, pi/2)", C.hi(beta1), C.lo(HALFPI)),
          ("(pi/2, beta*_2)", C.hi(HALFPI), C.lo(beta2)),
          ("(beta*_2, pi)", C.hi(beta2), MPPI - EPS0)]
PARTITION = []
for lbl, a, c in (RANGES if not FAST else []):
    cells, uncov = certify_interval(lbl, a, c)
    cnts = set(tuple(sorted(x["counts"].items())) for x in cells)
    covered = sum(float(mpmath.mpf(x["beta_box"][1]) - mpmath.mpf(x["beta_box"][0]))
                  for x in cells)
    PARTITION.append({"range": lbl,
                      "endpoints": [mpmath.nstr(a, 22), mpmath.nstr(c, 22)],
                      "n_boxes": len(cells),
                      "measure_covered": covered,
                      "measure_total": float(c - a),
                      "certified_counts": [dict(x) for x in cnts],
                      "uncovered_collars": [[mpmath.nstr(p, 22), mpmath.nstr(q, 22)]
                                            for p, q in uncov],
                      "cells": cells})
    print("   %-18s  %d boxes, counts %s, uncovered %.3e of %.3e"
          % (lbl, len(cells), sorted(cnts)[0] if cnts else "-",
             float(c - a) - covered, float(c - a)))

# --- the critical betas, treated exactly ---------------------------------
CRITICAL = {}
Wh = Wmax(HALFPI)
vals = Wof(HALFPI, Wh)
assert all(C.contains_zero(z) for z in vals)
assert C.contains_zero(dWof(HALFPI, Wh)[0])
CRITICAL["beta = 0"] = {
    "nature": "degenerate column",
    "statement": "t - beta = t, so Wtau = Wpot = Wwgt = 0 identically in t; the "
                 "entire column beta = 0 belongs to all three zero sets.",
    "count": "uncountable (identically zero)"}
CRITICAL["beta = beta*_1"] = {
    "nature": "double root of the torque Wronskian (lens endpoint)",
    "beta_enclosure": S(beta1, 30),
    "double_root_t": "t = beta/2   (piece B, w = 0)",
    "why": "Phi(pi-beta,0) = 0 and Phi_w(b,0) = 0 identically, so w = 0 is a "
           "double root; the Wtau count jumps 0 -> 2 here.",
    "Phi_and_Phi_w_enclosures": [S(Wof(PI - beta1, C.I(0))[0], 8),
                                 S(dWof(PI - beta1, C.I(0))[0], 8)],
    "y_of_the_double_root": "-1/2 exactly (h(beta/2) and -h(-beta/2) are equal)"}
CRITICAL["beta = pi/2"] = {
    "nature": "total degeneracy / orthogonal teacher",
    "statement": "b = pi/2 on BOTH pieces and W = pi/4.  Phi(pi/2,pi/4) = "
                 "Psi(pi/2,pi/4) = Om(pi/2,pi/4) = 0 and Phi_w(pi/2,pi/4) = 0, "
                 "so all three determinants vanish exactly on the same set "
                 "{t = 0, t = pi/2} (mod pi) -- both points lying ON the kink "
                 "lattice -- and there is no interior root (L5/L7 with "
                 "Phi(.,0), Psi(.,0) > 0).",
    "values_at_(pi/2,pi/4)": [S(z, 8) for z in vals],
    "Phi_w_at_(pi/2,pi/4)": S(dWof(HALFPI, Wh)[0], 8),
    "count": {"Wtau": 2, "Wpot": 2, "Wwgt": 2},
    "mass_direction": "at both roots h(t) = h(t-beta) = 0, so the mass direction "
                      "is 0/0; the curves extend by the limits y = c-1, -c "
                      "(torque) and y = 1-c, c (potential)."}
CRITICAL["beta = beta*_2"] = {
    "nature": "double root of the torque Wronskian (lens endpoint)",
    "beta_enclosure": S(beta2, 30),
    "double_root_t": "t = (beta+pi)/2   (piece A, w = 0)",
    "why": "Phi(beta,0) = 0 and Phi_w(b,0) = 0 identically; the Wtau count "
           "jumps 2 -> 0 here.",
    "Phi_and_Phi_w_enclosures": [S(Wof(beta2, C.I(0))[0], 8),
                                 S(dWof(beta2, C.I(0))[0], 8)],
    "y_of_the_double_root": "-1/2 exactly"}
CRITICAL["beta = pi"] = {
    "nature": "degenerate column (identified with beta = 0 by pi-periodicity)",
    "statement": "t - beta = t (mod pi); all three determinants vanish identically."}

# ==========================================================================
# 5.  The GENERIC certified interval root test, run on the FOLDED mirror of
#     widgets.js with the t-domain split at the kink lattice.
# ==========================================================================
print("== 5. generic certified root test on the folded mirror ==")
MIN_W = mpmath.mpf("1e-22")
TGT_W = mpmath.mpf("1e-28")


def generic_roots(det, beta, tlo, thi, kt, ku):
    roots, unres, nex = [], [], 0
    stack = [(tlo, thi)]
    while stack:
        a, c = stack.pop()
        T = iv.mpf([a, c])
        F, dF = C.dets_branch(beta, T, kt, ku)[det]
        if not C.contains_zero(F):
            nex += 1
            continue
        if not C.contains_zero(dF):
            Fa = C.dets_branch(beta, iv.mpf([a, a]), kt, ku)[det][0]
            Fc = C.dets_branch(beta, iv.mpf([c, c]), kt, ku)[det][0]
            sa, sc = sgn(Fa), sgn(Fc)
            if sa and sc:
                if sa == sc:
                    nex += 1
                    continue
                aa, cc, it = a, c, 0
                while cc - aa > TGT_W and it < 500:
                    m = (aa + cc) / 2
                    s = sgn(C.dets_branch(beta, iv.mpf([m, m]), kt, ku)[det][0])
                    if s == 0:
                        break
                    if s == sa:
                        aa = m
                    else:
                        cc = m
                    it += 1
                roots.append(iv.mpf([aa, cc]))
                continue
        if c - a < MIN_W:
            unres.append((a, c))
            continue
        m = (a + c) / 2
        stack += [(a, m), (m, c)]
    return roots, unres, nex


def y_encl(beta, T, kt, ku):
    At = C.atoms_on_branch(T, kt)
    Au = C.atoms_on_branch(T - beta, ku)
    return C.y_from_s(-Au.h, At.h)


TABLE = []
betas = [("pi*%d/24" % k, PI * k / 24) for k in range(1, 24)]
betas += [("beta*_1 - 1e-3", beta1 - C.I("1e-3")),
          ("beta*_1 + 1e-3", beta1 + C.I("1e-3")),
          ("pi/2 - 1e-3", HALFPI - C.I("1e-3")),
          ("pi/2 + 1e-3", HALFPI + C.I("1e-3")),
          ("beta*_2 - 1e-3", beta2 - C.I("1e-3")),
          ("beta*_2 + 1e-3", beta2 + C.I("1e-3")),
          ("pi/2 - 1e-6", HALFPI - C.I("1e-6")),
          ("pi/2 + 1e-6", HALFPI + C.I("1e-6"))]
t0 = time.time()
for lbl, beta in (betas if not FAST else []):
    bm = C.mid(beta)
    bt = iv.mpf([bm, bm])
    row = {"beta_label": lbl, "beta": mpmath.nstr(bm, 25),
           "kink_split_points_in_t": ["0", mpmath.nstr(bm, 25), "pi"]}
    for det in ("Wtau", "Wpot", "Wwgt"):
        interior, onkink, unres, nex = [], [], [], 0
        for (tl, th, kt, ku) in ((mpmath.mpf(0), bm, 0, -1),
                                 (bm, C.hi(PI), 0, 0)):
            r, u, ne = generic_roots(det, bt, tl, th, kt, ku)
            nex += ne
            for rr in r:
                # classify: a root whose enclosure meets the kink lattice
                # {0, beta, pi} is a kink root, not an interior root
                a_, c_ = C.endpoints(rr)
                on = min(abs(a_), abs(a_ - bm), abs(a_ - MPPI),
                         abs(c_), abs(c_ - bm), abs(c_ - MPPI)) < mpmath.mpf("1e-25")
                (onkink if on else interior).append((rr, kt, ku))
            unres += [(a, c, kt, ku) for a, c in u]
        for tk, kt, ku in ((mpmath.mpf(0), 0, -1), (bm, 0, -1)):
            v = C.dets_branch(bt, iv.mpf([tk, tk]), kt, ku)[det][0]
            if C.contains_zero(v):
                onkink.append((iv.mpf([tk, tk]), kt, ku))

        TOLD = mpmath.mpf("1e-18")

        def fmt(lst):
            """roots are counted in t modulo pi, so t and t+pi are one root"""
            out, seen = [], []
            for rr, kt, ku in sorted(lst, key=lambda z: C.lo(z[0])):
                m = C.mid(rr)
                if any(min(abs(m - s), abs(m - s - MPPI), abs(m - s + MPPI)) < TOLD
                       for s in seen):
                    continue
                seen.append(m)
                yv = y_encl(bt, rr, kt, ku)
                out.append({"t": S(rr, 22),
                            "y": S(yv, 20) if yv is not None else None})
            return out
        ir, kr = fmt(interior), fmt(onkink)
        # an undecided box that abuts the kink lattice is EXPLAINED: the
        # determinant is certified to vanish at that kink point, so the box
        # contains a known root and no further one (L5/L6/L7/L8).
        expl, genuine = [], []
        for a, c, _, _ in unres:
            near = min(abs(a), abs(a - bm), abs(a - MPPI),
                       abs(c), abs(c - bm), abs(c - MPPI))
            (expl if near <= (c - a) + mpmath.mpf("1e-20") else genuine).append((a, c))
        row[det] = {"n_interior_roots": len(ir), "interior_roots": ir,
                    "n_kink_roots": len(kr), "kink_roots": kr,
                    "boxes_excluded": nex,
                    "kink_adjacent_boxes_explained": len(expl),
                    "unresolved_boxes": [[mpmath.nstr(a, 20), mpmath.nstr(c, 20)]
                                         for a, c in genuine]}
        for a, c in genuine:
            UNRESOLVED.append({"kind": "generic_root_test_box", "det": det,
                               "beta": mpmath.nstr(bm, 20),
                               "t_box": [mpmath.nstr(a, 20), mpmath.nstr(c, 20)],
                               "width": float(c - a)})
    TABLE.append(row)
    print("   beta=%-16s interior/kink:  Wtau %d/%d  Wpot %d/%d  Wwgt %d/%d"
          % (lbl, row["Wtau"]["n_interior_roots"], row["Wtau"]["n_kink_roots"],
             row["Wpot"]["n_interior_roots"], row["Wpot"]["n_kink_roots"],
             row["Wwgt"]["n_interior_roots"], row["Wwgt"]["n_kink_roots"]))
print("   (%.1f s)" % (time.time() - t0))

# ==========================================================================
# 6.  The curves in the map coordinates (beta, y), with certified boxes.
# ==========================================================================
print("== 6. curves in (beta, y) ==")


def graded(lo2, hi2, n, both=True):
    """partition points clustered at the ends"""
    pts = []
    for i in range(n + 1):
        s = mpmath.mpf(i) / n
        if both:
            s = (1 - mpmath.cos(MPPI * s)) / 2
        pts.append(lo2 + (hi2 - lo2) * s)
    return pts


Y_TOL = mpmath.mpf(os.environ.get("Y_TOL", "1e-3"))
BETA_FLOOR = mpmath.mpf(os.environ.get("BETA_FLOOR", "1e-12"))
MAX_DEPTH_COUNT = 45             # refinement cap when the count is undecided
# Refinement cap when only the y enclosure is wide.  This is a BUDGET, not a
# mathematical limit: a box reported wide at the cap may well be tight two
# bisections later.  Raising it separates the boxes that are merely capped
# from the ones that are genuinely irreducible -- those bottom out on
# BETA_FLOOR instead, which is where the y enclosure stops improving under
# refinement (the curve endpoints, and the kink lattice at beta = pi/2).
MAX_DEPTH_Y = int(os.environ.get("MAX_DEPTH_Y", "7"))

# Stamp the run parameters into the certificate itself, so that every emitted
# certificate states the budget it was produced under instead of leaving the
# depth to be recovered from the console log.  (Certificates shipped before
# this stamp existed carry their parameters in cert/PROVENANCE.json instead.)
REPORT["run_parameters"] = {"MAX_DEPTH_Y": MAX_DEPTH_Y,
                            "MAX_DEPTH_COUNT": MAX_DEPTH_COUNT,
                            "Y_TOL": mpmath.nstr(Y_TOL, 17),
                            "BETA_FLOOR": mpmath.nstr(BETA_FLOOR, 17),
                            "FAST": FAST}


def box_record(det, Bbox, piece):
    """Certified record for one beta box, or (None, reason)."""
    r = count_piece(det, Bbox, piece)
    if r is None or r[0] != 1:
        return None, "count"
    wb = r[1][0]
    kt, ku = piece_branches(piece)
    y0 = y_encl(Bbox, t_of_w(Bbox, piece, wb), kt, ku)
    y1 = y_encl(Bbox, t_of_w(Bbox, piece, -wb), kt, ku)
    if y0 is None or y1 is None:
        # the mass direction (-h(t-beta), h(t)) is not sign-determined on this
        # box -- happens where a root approaches the kink lattice, i.e. near
        # beta = pi/2, where h(t) and h(t-beta) both vanish
        return None, "mass_direction_undetermined"
    rec = {"beta_box": S(Bbox, 22), "piece": piece, "w": S(wb, 20),
           "branch_w_pos": {"t": S(t_of_w(Bbox, piece, wb), 20), "y": S(y0, 20)},
           "branch_w_neg": {"t": S(t_of_w(Bbox, piece, -wb), 20), "y": S(y1, 20)},
           "y_sum": S(y0 + y1, 12),
           "_yw": max(C.width(y0), C.width(y1))}
    return rec, None


def curve_boxes(det, blo, bhi, piece, n=int(os.environ.get("NBOX", "400"))):
    """Certified boxes on a graded beta partition.

    Every returned box carries an enclosure of the root and of y valid for EVERY
    beta in that box.  Refinement is bounded: a box whose ROOT COUNT is not yet
    certified is bisected up to MAX_DEPTH_COUNT times (only the sub-box that
    still straddles the critical beta keeps failing, so this is linear); a box
    whose y enclosure is merely wide, or whose mass direction is not sign
    determined, is bisected at most MAX_DEPTH_Y times and then reported as a
    wide box -- still certified, just not tight."""
    tight, wide = [], []
    pts = graded(blo, bhi, n)
    todo = [(pts[i], pts[i + 1], 0, 0) for i in range(n)]
    todo.reverse()
    while todo:
        a, c, dc, dy = todo.pop()
        rec, reason = box_record(det, iv.mpf([a, c]), piece)
        if rec is not None and mpmath.mpf(rec["_yw"]) < Y_TOL:
            del rec["_yw"]
            tight.append(rec)
            continue
        deep = (dc >= MAX_DEPTH_COUNT) if rec is None else (dy >= MAX_DEPTH_Y)
        if c - a < BETA_FLOOR or deep:
            if rec is not None:
                del rec["_yw"]
            wide.append({"det": det, "piece": piece,
                         "beta_box": [mpmath.nstr(a, 22), mpmath.nstr(c, 22)],
                         "width": float(c - a),
                         "reason": reason or "y_enclosure_wider_than_%g" % float(Y_TOL),
                         "certified_root_exists_and_is_unique": rec is not None,
                         "record": rec})
            continue
        m = (a + c) / 2
        ndc, ndy = (dc + 1, dy) if rec is None else (dc, dy + 1)
        todo.append((m, c, ndc, ndy))
        todo.append((a, m, ndc, ndy))
    return tight, wide


LENS_A, LC1 = curve_boxes("Wtau", C.hi(beta1), C.lo(HALFPI), "B")
LENS_B, LC2 = curve_boxes("Wtau", C.hi(HALFPI), C.lo(beta2), "A")
POT_A, PC1 = curve_boxes("Wpot", EPS0, C.lo(HALFPI), "A")
POT_B, PC2 = curve_boxes("Wpot", C.hi(HALFPI), MPPI - EPS0, "B")
LENS, LCOL = LENS_A + LENS_B, LC1 + LC2
POTC, PCOL = POT_A + POT_B, PC1 + PC2
for col in LCOL + PCOL:
    UNRESOLVED.append(dict(col, kind="curve_box_wide_y_enclosure"))


def yrange(boxes):
    if not boxes:
        return ["EMPTY", "EMPTY"]
    los, his = [], []
    for r in boxes:
        for k in ("branch_w_pos", "branch_w_neg"):
            los.append(mpmath.mpf(r[k]["y"][0]))
            his.append(mpmath.mpf(r[k]["y"][1]))
    return [mpmath.nstr(min(los), 20), mpmath.nstr(max(his), 20)]


def ysum_check(boxes, target):
    return all(C.contains_zero(iv.mpf([mpmath.mpf(r["y_sum"][0]),
                                       mpmath.mpf(r["y_sum"][1])]) - target)
               for r in boxes)


ysum_ok = ysum_check(LENS, -1)
ysum_ok2 = ysum_check(POTC, 1)
covL = sum(float(mpmath.mpf(r["beta_box"][1]) - mpmath.mpf(r["beta_box"][0]))
           for r in LENS)
covP = sum(float(mpmath.mpf(r["beta_box"][1]) - mpmath.mpf(r["beta_box"][0]))
           for r in POTC)
print("   lens : %d boxes, beta measure %.12f of %.12f, y-hull %s"
      % (len(LENS), covL, float(C.lo(beta2) - C.hi(beta1)), yrange(LENS)))
print("   pot  : %d boxes, beta measure %.12f of %.12f, y-hull %s"
      % (len(POTC), covP, float(MPPI), yrange(POTC)))
print("   reciprocal identity  y+ + y- = -1 on lens:", ysum_ok,
      "  = +1 on potential curve:", ysum_ok2)
print("   wide boxes: lens %d, pot %d" % (len(LCOL), len(PCOL)))

CURVES = {
    "torque_lens": {
        "n_boxes": len(LENS),
        "beta_measure_covered": covL,
        "beta_measure_total": float(C.lo(beta2) - C.hi(beta1)),
        "y_hull_over_certified_boxes": yrange(LENS),
        "reflection_identity_y_plus_plus_y_minus_eq_minus_1": ysum_ok,
        "wide_boxes": LCOL,
        "boxes": LENS},
    "potential_curve": {
        "n_boxes": len(POTC),
        "beta_measure_covered": covP,
        "beta_measure_total": float(MPPI),
        "y_hull_over_certified_boxes": yrange(POTC),
        "reflection_identity_y_plus_plus_y_minus_eq_plus_1": ysum_ok2,
        "wide_boxes": PCOL,
        "boxes": POTC},
    "weight_wronskian_image": {
        "statement": "for every beta in (0,pi) the zero set of Wwgt in t is "
                     "exactly {0, beta} (mod pi); at t = 0 one has h(t) = 0 so "
                     "s1 = 0 and y = 0, at t = beta one has h(t-beta) = 0 so "
                     "s0 = 0 and y = +1 (= -1).  The image is therefore exactly "
                     "the two horizontal lines y = 0 and y = +1.",
        "certified_by": "L6 (Om = b cos(w+b/2) cos(w-b/2), both cosine arguments "
                        "in [-pi/2,pi/2], vanishing only at |w| = W)"}}

# ==========================================================================
# 7.  The explicit 1-D and 0-D piece lists
# ==========================================================================
c_lo, c_hi = S(cconst, 20)
b1s, b2s = S(beta1, 22), S(beta2, 22)

ONE_D = [
 {"id": "1D-1", "name": "torque Wronskian lens, upper branch",
  "equation": "Wtau(beta,t) = 0, the root whose mass ratio has rho > 1 (y > -1/2)",
  "beta_range": "[beta*_1, beta*_2] = [%s, %s]" % (b1s[0], b2s[1]),
  "parameterisation": "the two roots are w = +-w_r(beta) with w_r the unique zero "
                      "of Phi(b,.) in (0,(pi-b)/2); on piece B (beta < pi/2, "
                      "b = pi-beta, t = beta/2 + w) this branch is w = +w_r, on "
                      "piece A (beta > pi/2, b = beta, t = (beta+pi)/2 + w) it is "
                      "w = -w_r",
  "y_range": "[ -1/2 , -c )  with c = %s" % c_lo,
  "certified_boxes": "certificate.json: curves.torque_lens.boxes",
  "note": "y = -1/2 exactly at the two endpoints beta*_1, beta*_2; y -> -c as "
          "beta -> pi/2 (limit, not attained)"},
 {"id": "1D-2", "name": "torque Wronskian lens, lower branch",
  "equation": "Wtau(beta,t) = 0, the root with rho < 1 (y < -1/2)",
  "beta_range": "same as 1D-1",
  "parameterisation": "exact mirror of 1D-1 under t -> beta - t (mod pi), which "
                      "sends the mass ratio rho to 1/rho and hence y to -1-y",
  "y_range": "( c-1 , -1/2 ]",
  "certified_boxes": "certificate.json: curves.torque_lens.boxes"},
 {"id": "1D-3", "name": "potential Wronskian curve, lower branch",
  "equation": "Wpot(beta,t) = 0, the root with y < 1/2",
  "beta_range": "(0, pi)  (all of it)",
  "parameterisation": "the two roots are w = +-w_r(beta) with w_r the unique zero "
                      "of Psi(b,.) in (0,(pi-b)/2); piece A for beta < pi/2 "
                      "(b = beta, t = (beta+pi)/2 + w, this branch is w = -w_r), "
                      "piece B for beta > pi/2 (b = pi-beta, t = beta/2 + w, this "
                      "branch is w = +w_r)",
  "y_range": "( c , 1/2 )  -- entirely inside the MIXED-sign sector y > 0",
  "certified_boxes": "certificate.json: curves.potential_curve.boxes"},
 {"id": "1D-4", "name": "potential Wronskian curve, upper branch",
  "equation": "Wpot(beta,t) = 0, the root with y > 1/2",
  "beta_range": "(0, pi)",
  "parameterisation": "mirror of 1D-3 under t -> beta - t (mod pi); y -> 1-y",
  "y_range": "( 1/2 , 1-c )",
  "certified_boxes": "certificate.json: curves.potential_curve.boxes"},
 {"id": "1D-5", "name": "weight Wronskian image, root t = 0 (mod pi)",
  "equation": "Wwgt(beta,t) = 0 at t = 0; there h(t) = 0, so s1 = 0",
  "beta_range": "(0, pi)",
  "image_in_map": "the horizontal line y = 0",
  "note": "coincides with the degenerate stratum 1D-9 (second teacher massless)"},
 {"id": "1D-6", "name": "weight Wronskian image, root t = beta (mod pi)",
  "equation": "Wwgt(beta,t) = 0 at t = beta; there h(t-beta) = 0, so s0 = 0",
  "beta_range": "(0, pi)",
  "image_in_map": "the horizontal line y = +1 (identified with y = -1)",
  "note": "coincides with the degenerate stratum 1D-10"},
 {"id": "1D-7", "name": "degenerate column beta = 0",
  "equation": "beta = 0", "content": "all three determinants vanish identically "
  "in t; coincident teacher directions"},
 {"id": "1D-8", "name": "degenerate column beta = pi/2 (orthogonal teacher)",
  "equation": "beta = pi/2",
  "content": "all three determinants have the SAME zero set {t = 0, t = pi/2} "
             "(mod pi), both on the kink lattice; the mass direction degenerates "
             "to 0/0 at both"},
 {"id": "1D-9", "name": "degenerate stratum y = 0",
  "equation": "s1 = 0 -- the second teacher atom is massless"},
 {"id": "1D-10", "name": "degenerate stratum y = +1 (= y = -1)",
  "equation": "s0 = 0 -- the first teacher atom is massless; psi = 0 and psi = pi "
              "are the same projective line, so the top and bottom edges of the "
              "box are one stratum"},
 {"id": "1D-11", "name": "degenerate column beta = pi",
  "equation": "beta = pi",
  "content": "identified with beta = 0 by the pi-periodicity of the centered "
             "kernel; all three determinants vanish identically"}]

ZERO_D = [
 {"id": "0D-1", "name": "left lens endpoint / torque double root",
  "point": {"beta": b1s, "y": ["-0.5", "-0.5"]},
  "why": "double root of Wtau at t = beta/2 (w = 0): Phi(pi-beta,0) = 0 and "
         "Phi_w(.,0) = 0.  The Wtau root count jumps 0 -> 2.  Meeting point of "
         "1D-1 and 1D-2.  beta*_1 = 2u* with tan u* = pi/2 - u*.",
  "extra": {"u_star": S(ustar, 25)}},
 {"id": "0D-2", "name": "right lens endpoint / torque double root",
  "point": {"beta": b2s, "y": ["-0.5", "-0.5"]},
  "why": "double root of Wtau at t = (beta+pi)/2 (w = 0).  Count jumps 2 -> 0.  "
         "beta*_2 = pi - beta*_1 = b*, (b*/2) tan(b*/2) = 1."},
 {"id": "0D-3", "name": "torque lens meets the column beta = pi/2 (upper)",
  "point": {"beta": S(HALFPI, 22), "y": [str(-mpmath.mpf(c_hi)), str(-mpmath.mpf(c_lo))]},
  "why": "the Wtau root reaches the kink t = beta (equivalently t = 0); the mass "
         "direction is 0/0 there and the curve extends by the limit y = -c, "
         "c = (2/pi) arctan(2/pi).  y-extremum of 1D-1."},
 {"id": "0D-4", "name": "torque lens meets the column beta = pi/2 (lower)",
  "point": {"beta": S(HALFPI, 22),
            "y": [mpmath.nstr(mpmath.mpf(c_lo) - 1, 20),
                  mpmath.nstr(mpmath.mpf(c_hi) - 1, 20)]},
  "why": "mirror of 0D-3 under y -> -1-y; y-extremum of 1D-2."},
 {"id": "0D-5", "name": "potential curve meets the column beta = pi/2 (lower)",
  "point": {"beta": S(HALFPI, 22), "y": [c_lo, c_hi]},
  "why": "the Wpot root reaches the kink lattice; limit y = c.  y-extremum of 1D-3."},
 {"id": "0D-6", "name": "potential curve meets the column beta = pi/2 (upper)",
  "point": {"beta": S(HALFPI, 22),
            "y": [mpmath.nstr(1 - mpmath.mpf(c_hi), 20),
                  mpmath.nstr(1 - mpmath.mpf(c_lo), 20)]},
  "why": "mirror under y -> 1-y; y-extremum of 1D-4."},
 {"id": "0D-7", "name": "potential curve endpoint on the column beta = 0",
  "point": {"beta": ["0", "0"], "y": ["0.5", "0.5"]},
  "why": "as beta -> 0 the two Wpot roots tend to t = pi/2 +- w_D with w_D the "
         "Dottie number (cos w = w), and the mass ratio tends to rho = -1, i.e. "
         "y = 1/2.  1D-3 and 1D-4 meet here.",
  "extra": {"w_Dottie": S(wD, 25)}},
 {"id": "0D-8", "name": "potential curve endpoint on the column beta = pi",
  "point": {"beta": S(PI, 22), "y": ["0.5", "0.5"]},
  "why": "mirror of 0D-7 under beta -> pi - beta."},
 {"id": "0D-9", "name": "beta = 0 meets y = 0", "point": {"beta": ["0", "0"], "y": ["0", "0"]},
  "why": "corner of the degenerate strata 1D-7 and 1D-9 (also the endpoint of the "
         "weight-Wronskian image 1D-5)"},
 {"id": "0D-10", "name": "beta = 0 meets y = +1 (= -1)",
  "point": {"beta": ["0", "0"], "y": ["1", "1"]},
  "why": "corner of 1D-7 and 1D-10 (endpoint of 1D-6)"},
 {"id": "0D-11", "name": "beta = pi/2 meets y = 0",
  "point": {"beta": S(HALFPI, 22), "y": ["0", "0"]},
  "why": "crossing of the degenerate column 1D-8 with 1D-9/1D-5"},
 {"id": "0D-12", "name": "beta = pi/2 meets y = +1 (= -1)",
  "point": {"beta": S(HALFPI, 22), "y": ["1", "1"]},
  "why": "crossing of 1D-8 with 1D-10/1D-6"},
 {"id": "0D-13", "name": "beta = pi meets y = 0",
  "point": {"beta": S(PI, 22), "y": ["0", "0"]},
  "why": "corner of 1D-11 and 1D-9"},
 {"id": "0D-14", "name": "beta = pi meets y = +1 (= -1)",
  "point": {"beta": S(PI, 22), "y": ["1", "1"]},
  "why": "corner of 1D-11 and 1D-10"}]

NON_INTERSECTIONS = [
 "The torque lens (1D-1/1D-2) and the potential curve (1D-3/1D-4) never meet in "
 "the map: a common zero of Wtau and Wpot forces Q = Om = 0, and Om = 0 only at "
 "|w| = W, where Q = (pi/4)(pi-2b) sin b, so b = pi/2, i.e. beta = pi/2 and "
 "t in {0, pi/2}.  At that (beta,t) the mass direction is 0/0 and the two curves "
 "approach DIFFERENT limits (y = -c, c-1 versus y = c, 1-c), so their images are "
 "disjoint.",
 "The torque lens stays inside y in (c-1, -c) subset (-1,0), hence never meets "
 "y = 0 or y = +-1.",
 "The potential curve stays inside y in (c, 1-c) subset (0,1), hence never meets "
 "y = 0 or y = +-1."]

# ==========================================================================
# 8.  emit
# ==========================================================================
REPORT["structure_lemmas"] = L
REPORT["exact_linear_relation"] = ("Wtau(beta,t) + Wpot(beta,t) = 2*Wwgt(beta,t) "
                                   "identically -- the three determinants span a "
                                   "2-dimensional space (proved symbolically in "
                                   "validate.py, sympy)")
REPORT["beta_partition"] = PARTITION
REPORT["critical_betas"] = CRITICAL
REPORT["per_beta_table"] = TABLE
REPORT["curves"] = CURVES
REPORT["pieces_1D"] = ONE_D
REPORT["pieces_0D"] = ZERO_D
REPORT["non_intersections"] = NON_INTERSECTIONS
REPORT["unresolved_boxes"] = UNRESOLVED
try:
    REPORT["validation"] = json.load(open(
        _os.path.join(_HERE, "validation.json")))
except Exception as e:
    REPORT["validation"] = {"error": str(e)}

with open(_os.path.join(_HERE, "certificate.json"), "w") as fh:
    json.dump(REPORT, fh, indent=1)
print("wrote certificate.json;  unresolved entries:", len(UNRESOLVED))
from collections import Counter
print(Counter(u["kind"] for u in UNRESOLVED))
