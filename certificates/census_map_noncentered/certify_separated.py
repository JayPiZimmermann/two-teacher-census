"""
TASK 3 -- the SEPARATED stratum of the noncentered kernel, in TWO dimensions.

Object.  With  D = th0 - th1,  kap = pi,  and the balances of
LeanFormalization/Planar/GeneralTeachers/SeparatedBoundary.lean

    c0 = A(th1)/h(D),  c1 = -A(th0)/h(D),
    F0 = kap c0 + phi(D) c1 - P(th0),   F1 = phi(D) c0 + kap c1 - P(th1),
    separatedDet = F0(1,0) F1(0,1) - F0(0,1) F1(1,0),

the DIVISION-FREE determinant is  Dsep := h(D)^2 * separatedDet = a0 b1 - a1 b0
with the four column entries a0,b0,a1,b1 of noncentered.sep_det.  Dsep involves
only h and phi, which are C^1 everywhere (h' = slopeAtom, phi' = -h), so Dsep
is C^1 and BOTH its value and its gradient enclosures are sound across the kink
lattice; only the SECOND derivatives jump there.  Half-branch indices are
supplied explicitly whenever a box lies inside one half-branch of every one of
the five arguments th0, th1, th0-beta, th1-beta, D, and the multi-branch hull
is used otherwise.

CONDITIONING.  Dsep is ANTISYMMETRIC in (th0,th1) and vanishes to FOURTH order
on the diagonal D = 0 (mod 2pi): all four of a0,b0,a1,b1 vanish there.  The
test function used for EXCLUSION is therefore

    Psi := Dsep / h(D)^4          ( = separatedDet / h(D)^2 ),

which is O(1) as D -> 0 and has exactly the same zeros as separatedDet wherever
h(D) != 0.  h(D) = 0 exactly on D = 0 (mod pi).  Those two diagonals are
EXCLUDED BY MODELLING -- there the angular equations do not solve for the
masses, so separatedBalanceFst/Snd stop describing the separated stratum -- and
they are reported as their own stratum together with the collar half-width used.

The interval Newton test is run on Dsep itself (same zero set off h(D)=0),
using the mean-value enclosure of the value and the interval enclosure of the
partial derivative, so each box costs one gradient evaluation and one
thin-point evaluation.
"""
import math
import os

import mpmath
from mpmath import iv

import noncentered as J

PREC = 90
J.set_prec(PREC)
PIv = J.PI_IV()


def inter(x, y):
    a1, b1 = J.endpoints(x)
    a2, b2 = J.endpoints(y)
    a, b = max(a1, a2), min(b1, b2)
    if a > b:
        return x
    return iv.mpf([a, b])


def thin(x):
    return iv.mpf([x, x])


def branches_for(B, T0, T1):
    """Explicit half-branch indices of th0, th1, th0-B, th1-B, D, or None."""
    out = {}
    for key, X in (("n0", T0), ("n1", T1), ("n0b", T0 - B), ("n1b", T1 - B),
                   ("nD", T0 - T1)):
        ps = J.branch_pieces(X)
        if len(ps) != 1:
            return None
        out[key] = ps[0][1]
    return out


class BoxData(object):
    __slots__ = ("Dv", "D0", "D1", "Db", "Dc", "a0", "a1", "hD", "sAD",
                 "Psi", "ok")


def evaluate(B, T0, T1, br):
    """One gradient evaluation + one thin evaluation; returns BoxData or None
    when h(D) straddles 0 (then the box meets the h(D) = 0 stratum)."""
    Dv, D0, D1, Db, a0, a1, AD, sm, Dc = J.sep_det_grad_mv(B, T0, T1, br)
    hD, sAD = AD.h, AD.sA
    if J.contains_zero(hD):
        return None
    h2 = hD * hD
    h4 = h2 * h2
    r = BoxData()
    r.D0, r.D1, r.Db, r.a0, r.a1, r.hD, r.sAD = D0, D1, Db, a0, a1, hD, sAD
    r.Dv, r.Dc, r.Psi = Dv, Dc, Dv / h4
    return r


def newton(B, T0, T1, r, c):
    """Interval Newton on coordinate c using only the data in r.

    Returns ('unique', N) when Dsep has exactly one zero in the c-interval for
    EVERY beta in B and every value of the other angle in the box (a certified
    graph), 'empty' when it has none, or None."""
    X = T0 if c == 0 else T1
    a, b = J.endpoints(X)
    x0 = (a + b) / 2
    dP = r.D0 if c == 0 else r.D1
    if J.contains_zero(dP):
        return None, None
    b0_ = thin(J.mid(B))
    o0 = thin(J.mid(T1 if c == 0 else T0))
    O = T1 if c == 0 else T0
    dO = r.D1 if c == 0 else r.D0
    F = r.Dc + r.Db * (B - b0_) + dO * (O - o0)
    N = thin(x0) - F / dP
    na, nb = J.endpoints(N)
    if na > b or nb < a:
        return "empty", None
    if a < na and nb < b:
        return "unique", inter(N, X)
    return None, None


def refine(B, T0, T1, c, br, iters=30, tol=mpmath.mpf("1e-20")):
    X = T0 if c == 0 else T1
    for _ in range(iters):
        a, b = J.endpoints(X)
        if b - a < tol:
            break
        if c == 0:
            r = evaluate(B, X, T1, br)
        else:
            r = evaluate(B, T0, X, br)
        if r is None:
            break
        v, N = newton(B, X if c == 0 else T0, T1 if c == 0 else X, r, c)
        if v != "unique":
            break
        if J.width(N) >= J.width(X):
            break
        X = N
        if c == 0:
            T0 = X
        else:
            T1 = X
    return X


def halfdomain_bb(B, delta, minw, budget=400000, cell=0.4):
    """Cover the HALF domain (s, D) = (th1, th0-th1) with D in [delta, pi-delta],
    s in [0, 2pi], for a whole beta box B.

    That half domain carries the whole separated stratum: Dsep is ANTISYMMETRIC
    in (th0,th1), so the region D in (pi, 2pi) is the exact mirror image of
    D in (0, pi) under (th0,th1) -> (th1,th0), which also leaves the mass
    direction (-a1, a0) unchanged projectively.  The bands
    D in [0,delta] u [pi-delta, pi+delta] u [2pi-delta, 2pi] are the h(D) = 0
    stratum (plus a conditioning collar) and are reported separately."""
    two_pi = J.hi(2 * PIv)
    dpi = J.hi(PIv)
    scuts = {mpmath.mpf(0), two_pi}
    for base in (B, PIv, B + PIv, B - PIv, B + 2 * PIv):
        a, b = J.endpoints(base)
        for z in (a, b):
            if 0 < z < two_pi:
                scuts.add(z)
    scuts = sorted(scuts)
    sfine = []
    for i in range(len(scuts) - 1):
        a, b = scuts[i], scuts[i + 1]
        if b <= a:
            continue
        n = max(1, int(math.ceil(float(b - a) / cell)))
        for k in range(n):
            sfine.append((a + (b - a) * k / n, a + (b - a) * (k + 1) / n))
    dlo_, dhi_ = delta, dpi - delta
    nD = max(1, int(math.ceil(float(dhi_ - dlo_) / cell)))
    dfine = [(dlo_ + (dhi_ - dlo_) * k / nD, dlo_ + (dhi_ - dlo_) * (k + 1) / nD)
             for k in range(nD)]

    stack = [(a[0], a[1], b[0], b[1]) for a in sfine for b in dfine]
    domain_area = float(two_pi) * float(dhi_ - dlo_)
    excluded = mpmath.mpf(0)
    certified = mpmath.mpf(0)
    branch_boxes = []
    unresolved = []
    steps = 0
    n_eval = 0

    while stack:
        steps += 1
        if steps > budget:
            for bx in stack:
                unresolved.append({"box": [float(x) for x in bx],
                                   "area": float((bx[1] - bx[0]) * (bx[3] - bx[2])),
                                   "reason": "budget"})
            break
        s0, s1, d0, d1 = stack.pop()
        area = (s1 - s0) * (d1 - d0)
        T1 = iv.mpf([s0, s1])
        T0 = T1 + iv.mpf([d0, d1])
        br = branches_for(B, T0, T1)
        r = evaluate(B, T0, T1, br)
        n_eval += 1
        if r is None:
            unresolved.append({"box": [float(s0), float(s1), float(d0), float(d1)],
                               "area": float(area), "reason": "hD_straddles_zero"})
            continue
        if not J.contains_zero(r.Psi):
            excluded += area
            continue
        # gradients in (s, D):  d/ds = D0 + D1 ,  d/dD = D0
        gs = r.D0 + r.D1
        gd = r.D0
        got = None
        for c, g, X in ((0, gs, iv.mpf([s0, s1])), (1, gd, iv.mpf([d0, d1]))):
            if J.contains_zero(g):
                continue
            a, b = J.endpoints(X)
            x0 = (a + b) / 2
            og = gd if c == 0 else gs
            OX = iv.mpf([d0, d1]) if c == 0 else iv.mpf([s0, s1])
            o0 = thin(J.mid(OX))
            b0_ = thin(J.mid(B))
            F = r.Dc + r.Db * (B - b0_) + og * (OX - o0)
            N = thin(x0) - F / g
            na, nb = J.endpoints(N)
            if na > b or nb < a:
                got = "empty"
                break
            if a < na and nb < b:
                got = (c, inter(N, X))
                break
        if got == "empty":
            excluded += area
            continue
        if got is not None:
            c, N = got
            if c == 0:
                sN, dN = N, iv.mpf([d0, d1])
            else:
                sN, dN = iv.mpf([s0, s1]), N
            t1e = sN
            t0e = sN + dN
            rr = evaluate(B, t0e, t1e, branches_for(B, t0e, t1e))
            rec = {"box": [float(s0), float(s1), float(d0), float(d1)],
                   "graph_coord": "s" if c == 0 else "D", "area": float(area),
                   "s": [mpmath.nstr(J.lo(sN), 18), mpmath.nstr(J.hi(sN), 18)],
                   "D": [mpmath.nstr(J.lo(dN), 18), mpmath.nstr(J.hi(dN), 18)],
                   "th0": [mpmath.nstr(J.lo(t0e), 18), mpmath.nstr(J.hi(t0e), 18)],
                   "th1": [mpmath.nstr(J.lo(t1e), 18), mpmath.nstr(J.hi(t1e), 18)],
                   "branches": br}
            if rr is not None:
                y = J.y_from_s(-rr.a1, rr.a0)
                rec["y"] = None if y is None else [mpmath.nstr(J.lo(y), 18),
                                                   mpmath.nstr(J.hi(y), 18)]
                rec["y_width"] = None if y is None else J.width(y)
            else:
                rec["y"] = None
            branch_boxes.append(rec)
            certified += area
            continue
        if (s1 - s0) < minw and (d1 - d0) < minw:
            unresolved.append({"box": [float(s0), float(s1), float(d0), float(d1)],
                               "area": float(area), "reason": "minwidth"})
            continue
        if s1 - s0 >= d1 - d0:
            sm = (s0 + s1) / 2
            stack.append((s0, sm, d0, d1))
            stack.append((sm, s1, d0, d1))
        else:
            dm = (d0 + d1) / 2
            stack.append((s0, s1, d0, dm))
            stack.append((s0, s1, dm, d1))

    return {
        "beta": [mpmath.nstr(J.lo(B), 20), mpmath.nstr(J.hi(B), 20)],
        "precision_bits": PREC,
        "delta_collar": float(delta),
        "min_box_width": float(minw),
        "half_domain": "s in [0,2pi], D in [delta, pi-delta]",
        "half_domain_area": domain_area,
        "torus_area": float(two_pi) ** 2,
        "collar_area_hD_zero_stratum": float(two_pi) * 4 * float(delta),
        "excluded_area": float(excluded),
        "certified_branch_area": float(certified),
        "unresolved_area": float(sum(u.get("area", 0.0) for u in unresolved)),
        "closed_fraction_of_half_domain":
            float((excluded + certified) / mpmath.mpf(domain_area)),
        "n_branch_boxes": len(branch_boxes),
        "n_unresolved": len(unresolved),
        "n_evaluations": n_eval,
        "steps": steps,
        "branch_boxes": branch_boxes,
        "unresolved": unresolved[:400],
    }
