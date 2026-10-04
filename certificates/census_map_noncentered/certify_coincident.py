"""
TASK 2 -- the COINCIDENT stratum of the noncentered kernel.

Certified-complete enumeration of the roots in t (one 2pi period) of the
torque Wronskian Wtau and the potential Wronskian Wpot, for a certified
partition of beta in [0, 2pi).  Wwgt is determined by the exact identity
Wtau + Wpot = 2 Wwgt (validate.py (d)).

THE REDUCTION (see CERTIFICATE.md section 2).  With
   p = pi - (t mod 2pi) folded into [-pi,pi],  q = p - (x_t - x_u) ... i.e.
   p = pi - x_t,  q = pi - x_u  (x_t = t mod 2pi, x_u = (t-beta) mod 2pi),
   d = q - p  in  beta + 2pi Z,  m = (p+q)/2,
   e_t = sign(p), e_u = sign(q), E = e_t e_u:

   Wtau = |p q| * fhat,   Wpot = |p q| * ghat,   Wwgt = |p q| * S
   S(d,m)    = d * sinc(p) * sinc(q)              (sinc z = sin z / z, entire)
   fhat = S + sin d,   ghat = S - sin d,   sin d = sin beta.

|p q| >= 0 vanishes exactly at the kink lattice t in {0, beta, pi, beta+pi}.
The four kink points are the piece boundaries, and on each closed piece
fhat, ghat are ANALYTIC in m.  Values AT the piece boundaries:
   p = 0  (t = pi)      : fhat = 2 sin beta,  ghat = 0
   q = 0  (t = beta+pi) : fhat = 2 sin beta,  ghat = 0
   p = +-pi (t = 0)     : fhat = sin beta,    ghat = -sin beta
   q = +-pi (t = beta)  : fhat = sin beta,    ghat = -sin beta
so for sin beta != 0 the only boundary zeros are ghat at p=0 and q=0.
"""
import json
import math
import os
import sys

import mpmath
from mpmath import iv, mp

import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 200
J.set_prec(PREC)
MINW = mpmath.mpf("1e-11")          # minimum box width in m
MINW_BETA = mpmath.mpf("1e-13")     # minimum box width in beta

FACT = [mpmath.mpf(mpmath.factorial(k)) for k in range(0, 90)]


def _series(z, coeffs, M, tail):
    s = iv.mpf(0)
    pw = iv.mpf(1)
    z2 = z * z
    for c in coeffs:
        s = s + pw * c
        pw = pw * z2
    return s + iv.mpf([-tail, tail])


def sinc_iv(z):
    """sinc(z) = sin z / z, rigorous, removable singularity at 0."""
    a, b = J.endpoints(z)
    if a > 0 or b < 0:
        return iv.sin(z) / z
    M = float(max(abs(a), abs(b)))
    assert M <= 3.4, "sinc series branch needs |z| <= 3.4, got %s" % M
    K = 30
    coeffs = [((-1) ** k) / iv.mpf(str(FACT[2 * k + 1])) for k in range(K)]
    tail = float(mpmath.mpf(M) ** (2 * K) / FACT[2 * K + 1])
    return _series(z, coeffs, M, tail)


def sinc1_iv(z):
    """d/dz sinc(z) = (z cos z - sin z)/z^2."""
    a, b = J.endpoints(z)
    if a > 0 or b < 0:
        return (z * iv.cos(z) - iv.sin(z)) / (z * z)
    M = float(max(abs(a), abs(b)))
    assert M <= 3.4
    K = 30
    # sum_{k>=1} (-1)^k 2k z^{2k-1}/(2k+1)!  = z * sum_{j>=0} (-1)^(j+1) (2j+2) z^(2j)/(2j+3)!
    coeffs = [((-1) ** (j + 1)) * iv.mpf(2 * j + 2) / iv.mpf(str(FACT[2 * j + 3]))
              for j in range(K)]
    tail = float(mpmath.mpf(2 * K + 2) * mpmath.mpf(M) ** (2 * K) / FACT[2 * K + 3])
    return z * _series(z, coeffs, M, tail)


def sinc2_iv(z):
    """d2/dz2 sinc(z)."""
    a, b = J.endpoints(z)
    if a > 0 or b < 0:
        return ((z * z - 2) * iv.sin(z) + 2 * z * iv.cos(z)) / (z * z * z)
    M = float(max(abs(a), abs(b)))
    assert M <= 3.4
    K = 30
    # sum_{k>=1} (-1)^k 2k(2k-1) z^(2k-2)/(2k+1)!
    coeffs = [((-1) ** (j + 1)) * iv.mpf((2 * j + 2) * (2 * j + 1))
              / iv.mpf(str(FACT[2 * j + 3])) for j in range(K)]
    tail = float(mpmath.mpf((2 * K + 2) * (2 * K + 1)) * mpmath.mpf(M) ** (2 * K)
                 / FACT[2 * K + 3])
    return _series(z, coeffs, M, tail)


# Global rigorous bounds from sinc(z) = int_0^1 cos(zt) dt:
#   |sinc| <= 1, |sinc'| <= 1/2, |sinc''| <= 1/3.
_B0 = iv.mpf([-1, 1])
_B1 = iv.mpf([-mpmath.mpf(1) / 2, mpmath.mpf(1) / 2])
_B2 = iv.mpf([-mpmath.mpf(1) / 3, mpmath.mpf(1) / 3])


def inter(x, y):
    a1, b1 = J.endpoints(x)
    a2, b2 = J.endpoints(y)
    a, b = max(a1, a2), min(b1, b2)
    if a > b:
        return x
    return iv.mpf([a, b])


def sinc_all(z):
    """(sinc, sinc', sinc'') as a rigorous enclosure, intersected with the
    global bounds so that no cancellation in the direct forms can blow up."""
    return (inter(sinc_iv(z), _B0), inter(sinc1_iv(z), _B1),
            inter(sinc2_iv(z), _B2))


def _SS_raw(d, m):
    """(S, S_m, S_mm, S_d, S_md) with S = d sinc(p) sinc(q), p=m-d/2, q=m+d/2."""
    p = m - d / 2
    q = m + d / 2
    s0, s1, s2 = sinc_all(p)
    t0, t1, t2 = sinc_all(q)
    S = d * (s0 * t0)
    Sm = d * (s1 * t0 + s0 * t1)
    Smm = d * (s2 * t0 + 2 * (s1 * t1) + s0 * t2)
    Sd = s0 * t0 + d * ((s0 * t1 - s1 * t0) / 2)
    Smd = (s1 * t0 + s0 * t1) + (d / 2) * (s0 * t2 - s2 * t0)
    return S, Sm, Smm, Sd, Smd


def SS(d, m):
    """(S, S_m) with a two-dimensional mean-value refinement in (d, m)."""
    d0 = iv.mpf([J.mid(d), J.mid(d)])
    m0 = iv.mpf([J.mid(m), J.mid(m)])
    rd, rm = d - d0, m - m0
    S0, S0m, _, _, _ = _SS_raw(d0, m0)
    S, Sm, Smm, Sd, Smd = _SS_raw(d, m)
    Sm_mv = inter(Sm, S0m + Smd * rd + Smm * rm)
    S_mv = inter(S, S0 + Sd * rd + Sm_mv * rm)
    return S_mv, Sm_mv


def phi_pair(kind, d, sd, m):
    """(value, m-derivative) of fhat (kind='tau') or ghat (kind='pot')."""
    S, Sm = SS(d, m)
    return (S + sd, Sm) if kind == "tau" else (S - sd, Sm)


# --------------------------------------------------------------------------
# certified zero count of an analytic phi on a closed m-interval
# --------------------------------------------------------------------------

def _thin(x):
    return iv.mpf([x, x])


def count_zeros(kind, d, sd, innerA, innerB, sliverA, sliverB,
                zeroA, zeroB, minw=MINW, budget=1200):
    """Number of zeros of phi strictly inside the piece, for EVERY beta in the
    box, excluding the two known lattice zeros of ghat.

    The piece endpoints m_lo, m_hi are INTERVALS (beta is an interval), so the
    domain splits into
        sliverA = [lo(m_lo), hi(m_lo)]   -- may contain the true endpoint
        inner   = [hi(m_lo), lo(m_hi)]   -- inside the piece for every beta
        sliverB = [lo(m_hi), hi(m_hi)]
    On a sliver that carries a KNOWN zero (ghat at p=0 or q=0) we certify strict
    monotonicity, so the known zero is the only one there; on a sliver without a
    known zero we certify that phi does not vanish at all.
    Returns (count, unresolved, root_boxes)."""
    unresolved = []

    def val(lo_, hi_):
        return phi_pair(kind, d, sd, iv.mpf([lo_, hi_]))

    start, end = sliverA[0], sliverB[1]
    for (name, sl, known) in (("A", sliverA, zeroA), ("B", sliverB, zeroB)):
        a, b = sl
        if b <= a:
            if name == "A":
                start = a
            else:
                end = b
            continue
        if known:
            # peel a collar from the OUTER edge, at least as long as the sliver,
            # on which phi is certified strictly monotone; phi then has at most
            # one zero there, namely the known lattice zero inside the sliver.
            span = (innerB - innerA) + (b - a)
            ok = False
            for k in range(0, 60):
                delta = span / mpmath.mpf(2) ** k
                if delta < (b - a):
                    break
                box = (a, a + delta) if name == "A" else (b - delta, b)
                _, dv = val(box[0], box[1])
                if not J.contains_zero(dv):
                    if name == "A":
                        start = box[1]
                    else:
                        end = box[0]
                    ok = True
                    break
            if not ok:
                unresolved.append((float(a), float(b), "collar_%s_not_monotone" % name))
                return None, unresolved, []
        else:
            v, _ = val(a, b)
            if J.contains_zero(v):
                unresolved.append((float(a), float(b), "sliver_%s_not_excluded" % name))
                return None, unresolved, []
            if name == "A":
                start = b
            else:
                end = a

    lo_, hi_ = start, end
    if lo_ >= hi_:
        return 0, unresolved, []

    # ---- generic subdivision on the inner interval
    def safe_point(x, l, r):
        for k in range(40):
            for sgn in (0, 1, -1):
                y = x + sgn * (r - l) / mpmath.mpf(2) ** (k + 3)
                if not (l <= y <= r):
                    continue
                v, _ = val(y, y)
                if not J.contains_zero(v):
                    return y, v
        return None, None

    aa, va = safe_point(lo_, lo_, hi_)
    bb, vb = safe_point(hi_, lo_, hi_)
    if aa is None or bb is None or aa >= bb:
        unresolved.append((float(lo_), float(hi_), "endpoint_nonzero"))
        return None, unresolved, []
    for (x0, x1) in ((lo_, aa), (bb, hi_)):
        if x1 > x0:
            v, _ = val(x0, x1)
            if J.contains_zero(v):
                unresolved.append((float(x0), float(x1), "sliver_inner"))
                return None, unresolved, []

    total = 0
    root_boxes = []
    stack = [(aa, bb, va, vb)]
    steps = 0
    while stack:
        steps += 1
        if steps > budget:
            unresolved.append((float(stack[0][0]), float(stack[0][1]), "budget"))
            return None, unresolved, root_boxes
        x0, x1, v0, v1 = stack.pop()
        v, dv = val(x0, x1)
        if not J.contains_zero(v):
            continue
        if not J.contains_zero(dv):
            s0 = J.lo(v0) > 0
            s1 = J.lo(v1) > 0
            if s0 != s1:
                total += 1
                root_boxes.append((x0, x1))
            continue
        if x1 - x0 < minw:
            unresolved.append((float(x0), float(x1), "minwidth"))
            continue
        xm, vm = safe_point((x0 + x1) / 2, x0, x1)
        if xm is None:
            unresolved.append((float(x0), float(x1), "split"))
            continue
        stack.append((x0, xm, v0, vm))
        stack.append((xm, x1, vm, v1))
    return (total if not unresolved else None), unresolved, root_boxes


def refine_root(kind, d, sd, x0, x1, tol=mpmath.mpf("1e-25"), nmax=200):
    """Bisect a certified simple-root bracket [x0,x1] down to width tol."""
    v0 = phi_pair(kind, d, sd, iv.mpf([x0, x0]))[0]
    pos0 = J.lo(v0) > 0
    for _ in range(nmax):
        if x1 - x0 < tol:
            break
        xm = (x0 + x1) / 2
        vm = phi_pair(kind, d, sd, iv.mpf([xm, xm]))[0]
        if J.contains_zero(vm):
            break
        if (J.lo(vm) > 0) == pos0:
            x0 = xm
        else:
            x1 = xm
    return iv.mpf([x0, x1])


# --------------------------------------------------------------------------
# per-beta certification
# --------------------------------------------------------------------------

def piece_endpoint_kind(pc):
    """For each end of the m-interval, whether it is a p=0 / q=0 point."""
    out = []
    for m in (pc["m_lo"], pc["m_hi"]):
        p = m - pc["d"] / 2
        q = m + pc["d"] / 2
        out.append(J.contains_zero(p) or J.contains_zero(q))
    return out


def certify_beta_box(B, minw=MINW):
    """B: interval of beta (must avoid 0, pi, 2pi).  Returns dict."""
    sd = iv.sin(B)
    if J.contains_zero(sd):
        return {"ok": False, "reason": "sin beta straddles 0"}
    pcs = J.pieces_of(B)
    res = {"pieces": [], "Wtau_interior": 0, "Wpot_interior": 0,
           "unresolved": []}
    for idx, pc in enumerate(pcs):
        d = pc["d"]
        A = J.hi(pc["m_lo"])       # inner enclosure endpoints, so that the
        Bx = J.lo(pc["m_hi"])      # scanned interval is inside every beta's piece
        outerA = J.lo(pc["m_lo"])
        outerB = J.hi(pc["m_hi"])
        za, zb = piece_endpoint_kind(pc)
        rec = {"idx": idx, "t_lo": float(J.mid(pc["t_lo"])),
               "t_hi": float(J.mid(pc["t_hi"])), "d": float(J.mid(d)),
               "E": pc["E"], "m_lo": float(J.mid(pc["m_lo"])),
               "m_hi": float(J.mid(pc["m_hi"])), "zero_ends": [za, zb]}
        slA = (J.lo(pc["m_lo"]), J.hi(pc["m_lo"]))
        slB = (J.lo(pc["m_hi"]), J.hi(pc["m_hi"]))
        for kind, key in (("tau", "Wtau"), ("pot", "Wpot")):
            zA = za and kind == "pot"
            zB = zb and kind == "pot"
            c, un, rb = count_zeros(kind, d, sd, A, Bx, slA, slB, zA, zB, minw)
            rec[key] = c
            rec[key + "_unresolved"] = un
            rec[key + "_root_boxes"] = rb
            if c is None:
                res["unresolved"].append({"piece": idx, "det": key, "boxes": un})
            else:
                res[key + "_interior"] += c
        res["pieces"].append(rec)
    res["ok"] = not res["unresolved"]
    # total root counts in t over [0,2pi): the two lattice roots t=pi, t=beta+pi
    # are roots of BOTH Wtau (touching, |p|-kink, no sign change) and Wpot
    # (sign change with vanishing t-derivative)
    if res["ok"]:
        res["Wtau_total"] = 2 + res["Wtau_interior"]
        res["Wpot_total"] = 2 + res["Wpot_interior"]
    return res


def lens_y(B, pc, mbox):
    """(s0,s1) and y enclosure of a root located in mbox on piece pc."""
    d = pc["d"]
    p = mbox - d / 2
    q = mbox + d / 2
    et = 1 if pc["nt"] % 2 == 0 else -1
    eu = 1 if pc["nu"] % 2 == 0 else -1
    ht = (et * p) * iv.sin(p)
    hu = (eu * q) * iv.sin(q)
    return -hu, ht, J.y_from_s(-hu, ht)


# --------------------------------------------------------------------------
# beta* : the double root of Wtau, tan(u) = pi - u with u = beta/2
# --------------------------------------------------------------------------

def certify_beta_star():
    """g(u) = sin u - (pi-u) cos u has a unique zero in (0, pi/2);
    beta*_1 = 2u*, beta*_2 = 2pi - 2u*."""
    def g(u):
        return iv.sin(u) - (J.PI_IV() - u) * iv.cos(u)

    def gp(u):
        return 2 * iv.cos(u) + (J.PI_IV() - u) * iv.sin(u)

    lo_, hi_ = mpmath.mpf("1.11"), mpmath.mpf("1.12")
    br = iv.mpf([lo_, hi_])
    assert J.lo(gp(br)) > 0, "g not certified increasing on the bracket"
    assert J.hi(g(iv.mpf([lo_, lo_]))) < 0 < J.lo(g(iv.mpf([hi_, hi_])))
    for _ in range(400):
        m = (lo_ + hi_) / 2
        v = g(iv.mpf([m, m]))
        if J.contains_zero(v):
            break
        if J.hi(v) < 0:
            lo_ = m
        else:
            hi_ = m
        if hi_ - lo_ < mpmath.mpf("1e-55"):
            break
    u = iv.mpf([lo_, hi_])
    return u, 2 * u, 2 * J.PI_IV() - 2 * u


if __name__ == "__main__":
    out = {"precision_bits": PREC, "min_box_width_m": float(MINW)}
    u, b1, b2 = certify_beta_star()
    out["beta_star"] = {
        "u_star": [mpmath.nstr(J.lo(u), 30), mpmath.nstr(J.hi(u), 30)],
        "characterisation": "tan u = pi - u, u in (0, pi/2); equivalently "
                            "fhat(d,0) = 0 on the piece t in [0,beta] with "
                            "d = beta - 2pi",
        "beta_star_1": [mpmath.nstr(J.lo(b1), 30), mpmath.nstr(J.hi(b1), 30)],
        "beta_star_2": [mpmath.nstr(J.lo(b2), 30), mpmath.nstr(J.hi(b2), 30)],
        "sum_minus_2pi": [mpmath.nstr(x, 8) for x in
                          J.endpoints(b1 + b2 - 2 * J.PI_IV())],
    }
    print("beta*_1 in [%s, %s]" % (mpmath.nstr(J.lo(b1), 25), mpmath.nstr(J.hi(b1), 25)))
    print("beta*_2 in [%s, %s]" % (mpmath.nstr(J.lo(b2), 25), mpmath.nstr(J.hi(b2), 25)))
    print(json.dumps(out["beta_star"], indent=1))
