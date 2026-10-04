"""
Certified core for the arrangement lane.

Everything here goes through mpmath.iv at 200 bits (directed rounding), reusing
the kernel module of the certificate lane (cert/centered.py).  float64 is never
used on a reported path.

The object certified here is the SIGN CHART of a teacher (beta, y):

    * the complete list of torque roots  t in [0, pi)  of
          F(t) = s0 h(t) + s1 h(t - beta),      (s0,s1) = (sin psi, cos psi),
          psi = (y+1) pi / 2
      (h = couplingH; F = -Trq of widgets.js), and
    * at each root, the certified signs of the three mass-free determinants
          Wtau, Wpot, Wwgt   of certificate section 0.

Why that is the census.  At a torque root the mass direction is
(s0,s1) = kappa * (-h(t-beta), h(t)) for a nonzero real kappa, hence

    P(t)   = s0 phi(t)   + s1 phi(t-beta)   = -kappa * Wpot(beta,t)
    W(t)   = s0|sin t|   + s1|sin(t-beta)|  = -kappa * Wwgt(beta,t)
    tau(t) = P(t) - 2 W(t)                  = +kappa * Wtau(beta,t)

(the last using Wtau + Wpot = 2 Wwgt).  widgets.js `coincidenceTypeAt` calls a
coincidence row a spurious minimum iff  tau*P > 0  and  W*a*b*tau < 0, so with
a*b > 0 for the same-sign split and a*b < 0 for the opposite split:

    trap present   <=>  tau*P > 0      <=>  Wtau * Wpot < 0     (kappa cancels)
    label @positive <=> W*tau < 0      <=>  Wtau * Wwgt > 0
    label @mixed    <=> W*tau > 0      <=>  Wtau * Wwgt < 0

so the census signature is a function of the sign chart alone.
"""
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "cert"))

import mpmath
from mpmath import iv, mp
import centered as C

C.set_prec(200)

PI = C.PI_IV
SPLIT = mpmath.mpf(0.5173648177666931)


def ivx(a, b=None):
    return C.I(a, b)


def sgn(x):
    """Certified strict sign of an interval, or 0 if it straddles."""
    a, b = C.endpoints(x)
    if a > 0:
        return 1
    if b < 0:
        return -1
    return 0


def masses(y):
    """(s0, s1) enclosure for the map coordinate y (exact input)."""
    psi = (C.I(y) + 1) * PI() / 2
    return iv.sin(psi), iv.cos(psi)


class Root(object):
    __slots__ = ("piece", "kt", "ku", "T", "wtau", "wpot", "wwgt")


def _F(B, s0, s1, T, kt, ku):
    At = C.atoms_on_branch(T, kt)
    Au = C.atoms_on_branch(T - B, ku)
    return s0 * At.h + s1 * Au.h, s0 * At.sA + s1 * Au.sA


def roots_on_piece(B, s0, s1, tlo, thi, kt, ku, tol=None, maxdepth=90):
    """Certified-complete isolation of the zeros of F on the CLOSED piece
    [tlo, thi] (given as mp numbers), which must lie inside branch kt for t and
    branch ku for t - B.

    Returns (roots, undecided) where roots is a list of tight enclosures and
    undecided is a list of boxes the test could not resolve (empty = complete).
    """
    if tol is None:
        tol = mpmath.mpf(2) ** -120
    todo = [(tlo, thi, 0)]
    roots, undecided = [], []
    while todo:
        a, b, d = todo.pop()
        T = iv.mpf([a, b])
        F, dF = _F(B, s0, s1, T, kt, ku)
        if sgn(F) != 0:
            continue                                  # certified no root
        if b - a <= tol:
            roots.append(T)                           # tight enough
            continue
        if sgn(dF) != 0:
            Fa, _ = _F(B, s0, s1, iv.mpf([a, a]), kt, ku)
            Fb, _ = _F(B, s0, s1, iv.mpf([b, b]), kt, ku)
            sa, sb = sgn(Fa), sgn(Fb)
            if sa != 0 and sb != 0:
                if sa == sb:
                    continue                          # monotone, no sign change
                # strictly monotone with a sign change: exactly one root.
                # bisect, then finish with interval Newton (which also handles a
                # root that lands exactly on a bisection midpoint).
                lo, hi = a, b
                for _ in range(60):
                    if hi - lo <= tol:
                        break
                    m = (lo + hi) / 2
                    Fm, _ = _F(B, s0, s1, iv.mpf([m, m]), kt, ku)
                    sm = sgn(Fm)
                    if sm == 0:
                        break
                    if sm == sa:
                        lo = m
                    else:
                        hi = m
                Tr = iv.mpf([lo, hi])
                for _ in range(20):
                    m = C.mid(Tr)
                    Fm, _ = _F(B, s0, s1, iv.mpf([m, m]), kt, ku)
                    _, dT = _F(B, s0, s1, Tr, kt, ku)
                    if sgn(dT) == 0:
                        break
                    N = C.I(m) - Fm / dT
                    nl, nh = C.endpoints(N)
                    tl, th = C.endpoints(Tr)
                    nl, nh = max(nl, tl), min(nh, th)
                    if nl > nh:
                        break
                    if nh - nl >= th - tl:
                        break
                    Tr = iv.mpf([nl, nh])
                    if nh - nl <= tol:
                        break
                roots.append(Tr)
                continue
        if d >= maxdepth:
            undecided.append((a, b))
            continue
        # split off-centre: exact roots of the symmetric strata (y = -1/2,
        # t = beta/2, t = beta/2 + pi/2, ...) land exactly on midpoints, and a
        # root sitting on a subdivision point can never be bracketed.
        m = a + (b - a) * SPLIT
        todo.append((a, m, d + 1))
        todo.append((m, b, d + 1))
    roots.sort(key=lambda r: C.endpoints(r)[0])
    # merge enclosures that touch (can only happen at a subdivision point)
    merged = []
    for r in roots:
        if merged and C.endpoints(merged[-1])[1] >= C.endpoints(r)[0]:
            merged[-1] = C.hull(merged[-1], r)
        else:
            merged.append(r)
    return merged, undecided


def sign_chart(beta, y, tol=None):
    """Certified sign chart at the exact map point (beta, y).

    beta, y are exact inputs (str / int / mpf).  Returns a dict.
    """
    B = C.I(beta)
    s0, s1 = masses(y)
    bl, bh = C.endpoints(B)
    pil, pih = C.endpoints(PI())
    assert bl > 0 and bh < pil, "beta must be interior"
    assert bl == bh, "beta must be exact"
    b = bl
    pi_lo = pil

    out = {"beta": str(b), "y": str(y), "roots": [], "undecided": []}
    # piece B : t in [0, beta],  branch kt = 0, ku = -1
    # piece A : t in [beta, pi], branch kt = 0, ku = 0
    for name, (a0, a1, kt, ku) in (
            ("B", (mpmath.mpf(0), b, 0, -1)),
            ("A", (b, pi_lo, 0, 0))):
        rs, und = roots_on_piece(B, s0, s1, a0, a1, kt, ku, tol=tol)
        out["undecided"] += [(name, str(x), str(z)) for x, z in und]
        for T in rs:
            D = C.dets_branch(B, T, kt, ku)
            r = {"piece": name,
                 "t": [str(C.endpoints(T)[0]), str(C.endpoints(T)[1])],
                 "t_mid": mpmath.nstr(C.mid(T), 25),
                 "sgn_Wtau": sgn(D["Wtau"][0]),
                 "sgn_Wpot": sgn(D["Wpot"][0]),
                 "sgn_Wwgt": sgn(D["Wwgt"][0]),
                 "Wtau": mpmath.nstr(C.mid(D["Wtau"][0]), 12),
                 "Wpot": mpmath.nstr(C.mid(D["Wpot"][0]), 12),
                 "Wwgt": mpmath.nstr(C.mid(D["Wwgt"][0]), 12)}
            out["roots"].append(r)
    out["nroots"] = len(out["roots"])
    out["census"] = census_of(out["roots"])
    out["certified"] = (not out["undecided"]) and all(
        r["sgn_Wtau"] and r["sgn_Wpot"] and r["sgn_Wwgt"] for r in out["roots"])
    return out


def census_of(roots):
    """widgets.js censusSignature, reconstructed from the certified sign chart."""
    items = ["fit:global"]
    for r in roots:
        st, sp, sw = r["sgn_Wtau"], r["sgn_Wpot"], r["sgn_Wwgt"]
        if st == 0 or sp == 0 or sw == 0:
            return None                        # not certified
        if st * sp < 0:                        # trap present
            items.append("coincident:trap@" + ("positive" if st * sw > 0 else "mixed"))
    return " | ".join(sorted(items))


# --------------------------------------------------------------------------
# certified branch values of the two curves at a given beta
# --------------------------------------------------------------------------

def _redu_pieces(beta):
    """(piece, b, sigma, W, t_of_w) for the two smooth pieces at exact beta."""
    B = C.I(beta)
    return [("A", B, 1, (PI() - B) / 2, lambda w: (B + PI()) / 2 + w),
            ("B", PI() - B, -1, (PI() - (PI() - B)) / 2, lambda w: B / 2 + w)]


def _bisect_w(f, lo, hi, tol):
    """Unique root of f (certified strictly monotone on [lo,hi]) by bisection."""
    flo, fhi = f(iv.mpf([lo, lo])), f(iv.mpf([hi, hi]))
    slo, shi = sgn(flo), sgn(fhi)
    assert slo != 0 and shi != 0 and slo != shi, "no certified bracket"
    while hi - lo > tol:
        m = (lo + hi) / 2
        sm = sgn(f(iv.mpf([m, m])))
        if sm == 0:
            break
        if sm == slo:
            lo = m
        else:
            hi = m
    return iv.mpf([lo, hi])


def curve_branches(beta, which, tol=None):
    """Certified y-enclosures of the two branches of the torque lens
    (which='torque') or of the potential curve (which='potential') at exact beta.

    Returns a list of dicts, or [] if the curve is absent at that beta.
    """
    if tol is None:
        tol = mpmath.mpf(2) ** -120
    B = C.I(beta)
    half = PI() / 2
    out = []
    for piece, bI, sigma, WI, t_of_w in _redu_pieces(beta):
        blo, bhi = C.endpoints(bI)
        # which piece carries the roots:
        #   Psi has its two zeros on the piece with b < pi/2
        #   Phi has its two zeros on the piece with pi/2 < b < b*
        if which == "potential":
            if not (bhi < C.endpoints(half)[0]):
                continue
            f = lambda w: C.redu(bI, w)[1]
            lo = C.endpoints(bI / 2)[1]          # certified lower bracket b/2
        else:
            if not (blo > C.endpoints(half)[1]):
                continue
            bstar_hi = mpmath.mpf("1.7206671780387595249678")
            if not (bhi < bstar_hi):
                continue
            f = lambda w: C.redu(bI, w)[0]
            lo = mpmath.mpf(0)
        hi = C.endpoints(WI)[0]                  # certified inside the piece
        Wroot = _bisect_w(f, lo, hi, tol)
        for s in (1, -1):
            w = Wroot if s == 1 else -Wroot
            T = t_of_w(w)
            s0v, s1v = C.y_of_root(B, T)
            yv = C.y_from_s(s0v, s1v)
            out.append({"piece": piece, "sign": s,
                        "w": mpmath.nstr(C.mid(Wroot), 25),
                        "t": mpmath.nstr(C.mid(T), 25),
                        "y_lo": (mpmath.nstr(C.endpoints(yv)[0], 25) if yv is not None else None),
                        "y_hi": (mpmath.nstr(C.endpoints(yv)[1], 25) if yv is not None else None),
                        "y_width": (C.width(yv) if yv is not None else None)})
    return out
