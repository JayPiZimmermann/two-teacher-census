"""The 2-D mass-free separated sweep in the WELL-CONDITIONED direction.

Same object as certify_separated.py -- a certified covering of the half domain
(s, D) = (th1, th0-th1) by boxes that either certifiably contain NO zero of the
mass-free separated determinant, or certifiably carry a unique branch of its
zero curve -- but with the RADIAL-ROW elimination (noncentered.wc_columns /
wc_det_grad_mv):

    DsepWC = a0 b1 - a1 b0,   a,b the massDet-cleared torque residuals at the
                              two unit teachers,

whose clearing determinant massDet = pi^2 - phiJ(D)^2 vanishes ONLY at D = 0
(mod 2pi).  Two consequences over the h(D)-cleared form of certify_separated.py:

  * the ANTIPODAL stratum D = pi is inside the certified domain (the old form
    degenerates there: h(pi) = 0), so the sweep runs on D in [delta, pi], not
    [delta, pi - delta];
  * the coincidence collar delta can be taken small, since the elimination is
    regular for every D != 0.

A separated critical configuration of ANY teacher (s0,s1) != (0,0) with
massDet != 0 satisfies DsepWC = 0 (linearity of the cleared residuals in the
teacher masses -- the same argument as separatedDet_eq_zero_of_
isSeparatedBalanced, run on the radial-solved residual), so a box excluded by
DsepWC != 0 contains no separated critical configuration of any teacher.
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
    out = {}
    for key, X in (("n0", T0), ("n1", T1), ("n0b", T0 - B), ("n1b", T1 - B),
                   ("nD", T0 - T1)):
        ps = J.branch_pieces(X)
        if len(ps) != 1:
            return None
        out[key] = ps[0][1]
    return out


class BoxData(object):
    __slots__ = ("Dv", "D0", "D1", "Db", "Dc", "a0", "a1", "M")


def evaluate(B, T0, T1, br):
    """Exclusion/Newton data from the REDUCED determinant Dred (noncentered.
    dred_grad_mv).  DsepWC = massDet * Dred (sympy-exact), and massDet > 0 for
    D != 0 (mod 2pi), so {Dred = 0} is the same zero set on the sweep domain
    and its enclosures avoid the massDet clearing factor's dependency."""
    Dv, D0, D1, Db, AD, M, sm, Dc = J.dred_grad_mv(B, T0, T1, br)
    if J.contains_zero(M):
        return None
    r = BoxData()
    r.Dv, r.D0, r.D1, r.Db, r.Dc, r.M = Dv, D0, D1, Db, Dc, M
    r.a0 = r.a1 = None
    return r


def kernel_vector(B, T0, T1, br):
    """(a0, a1) of the WC columns, for the mass direction of a branch box."""
    C, AD, M, sm = J.wc_columns(B, T0, T1, br)
    return C["a0"], C["a1"]


def halfdomain_bb(B, delta, minw, budget=400000, cell=0.4, d_hi=None):
    """Cover the half domain (s, D), s in [0, 2pi], D in [delta, pi], for the
    whole beta box B.  Returns the same record shape as
    certify_separated.halfdomain_bb."""
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
    dlo_ = delta
    dhi_ = dpi if d_hi is None else d_hi
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
        Dm = iv.mpf([d0, d1])
        T0 = T1 + Dm
        br = branches_for(B, T0, T1)
        r = evaluate(B, T0, T1, br)
        n_eval += 1
        if r is None:
            unresolved.append({"box": [float(s0), float(s1), float(d0), float(d1)],
                               "area": float(area), "reason": "massDet_straddles_zero"})
            continue
        if not J.contains_zero(r.Dv):
            excluded += area
            continue
        # gradients in (s, D):  d/ds = D0 + D1 ,  d/dD = D0
        gs = r.D0 + r.D1
        gd = r.D0
        got = None
        for c, g in ((0, gs), (1, gd)):
            if J.contains_zero(g):
                continue
            X = iv.mpf([s0, s1]) if c == 0 else iv.mpf([d0, d1])
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
                   "D": [mpmath.nstr(J.lo(dN), 18), mpmath.nstr(J.hi(dN), 18)]}
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
            m = (s0 + s1) / 2
            stack.append((s0, m, d0, d1))
            stack.append((m, s1, d0, d1))
        else:
            m = (d0 + d1) / 2
            stack.append((s0, s1, d0, m))
            stack.append((s0, s1, m, d1))

    return {
        "beta": [mpmath.nstr(J.lo(B), 20), mpmath.nstr(J.hi(B), 20)],
        "direction": "well-conditioned (radial rows / massDet)",
        "precision_bits": PREC,
        "delta_collar": float(delta),
        "min_box_width": float(minw),
        "half_domain": "s in [0,2pi], D in [delta, pi]",
        "half_domain_area": domain_area,
        "collar_area_massDet_zero_stratum": float(two_pi) * float(delta),
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


if __name__ == "__main__":
    import json
    import sys
    import time
    bc = mpmath.mpf(sys.argv[1])
    hw = mpmath.mpf(sys.argv[2]) if len(sys.argv) > 2 else mpmath.mpf(0)
    delta = mpmath.mpf(sys.argv[3]) if len(sys.argv) > 3 else mpmath.mpf("0.05")
    minw = mpmath.mpf(sys.argv[4]) if len(sys.argv) > 4 else mpmath.mpf("2e-4")
    budget = int(sys.argv[5]) if len(sys.argv) > 5 else 400000
    B = iv.mpf([bc - hw, bc + hw])
    t0 = time.time()
    r = halfdomain_bb(B, delta, minw, budget=budget)
    r["seconds"] = time.time() - t0
    print("beta %s +-%s: closed %.6f, %d branch boxes, %d unresolved "
          "(area %.3e), %d evals, %.0fs"
          % (sys.argv[1], mpmath.nstr(hw, 5), r["closed_fraction_of_half_domain"],
             r["n_branch_boxes"], r["n_unresolved"], r["unresolved_area"],
             r["n_evaluations"], r["seconds"]), flush=True)
    tag = sys.argv[6] if len(sys.argv) > 6 else sys.argv[1]
    with open("separated_wc_%s.json" % tag, "w") as fh:
        json.dump(r, fh)
