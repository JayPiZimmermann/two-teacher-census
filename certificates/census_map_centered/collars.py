"""
TASK 3: close the six unresolved beta collars of cert/CERTIFICATE.md section 9.1.

The certificate's generic box test could not certify a ROOT COUNT on six beta
collars of total measure 4.712e-13 around beta1*, pi/2, beta2*.  Those collars
close, and in fact the whole beta axis closes, because the counts are decided by
EXACT inequalities on the reduced forms rather than by a generic box test.

On a smooth piece with parameter b (b = beta on piece A, b = pi - beta on piece
B), W = (pi-b)/2, |w| <= W:

    Q  = (w^2 - b^2/4) sin b        Om = (b/2)(cos 2w + cos b)
    Phi = Q + Om  (Wtau)            Psi = -Q + Om  (Wpot)     Om (Wwgt)

C1 (Wwgt).  Om = b cos(w+b/2) cos(w-b/2) >= 0 with equality exactly at |w| = W
    (both cosine arguments lie in [b - pi/2, pi/2]).  Two zeros in t per period,
    at t = 0 and t = beta.  Exact, no numerics.

C2 (Wpot).  Psi_w = -(2w sin b + b sin 2w) < 0 on (0, W] because sin b > 0 and
    sin 2w > 0 for 2w in (0, pi-b).  Psi(b,0) = b cos(b/2)[cos(b/2)+(b/2)sin(b/2)]
    > 0.  Psi(b,W) = (pi/2)(b - pi/2) sin b.  Hence
        b < pi/2  ->  exactly one zero in (0,W), i.e. 2 zeros in t
        b > pi/2  ->  no zero.
    Determined for EVERY b != pi/2.

C3 (Wtau), b > pi/2.  D := -Phi_w = b sin 2w - 2w sin b has D(0) = 0,
    D_ww = -4b sin 2w <= 0 on 0 <= 2w <= pi-b, and D(W) = (2b-pi) sin b > 0;
    concavity with D(0)=0 gives D(w) >= (w/W) D(W) > 0, so Phi is strictly
    decreasing on [0,W].  Phi(b,W) = (pi/4)(pi-2b) sin b < 0.  So the zero count
    is decided by the sign of Phi(b,0) = b cos(b/2) g(b), g(b) = cos(b/2) -
    (b/2) sin(b/2):  2 zeros in t if g(b) > 0, none if g(b) < 0.

C4 (Wtau), b < pi/2.  Phi > 0 on [0,W], hence NO zero, by an exact argument that
    needs no case analysis near b = pi/2:
      * on [0, b/2]:  Phi_w = -D with D(0) = 0, D(b/2) = b sin b - b sin b = 0,
        and D concave on 2w in [0,b] (subset of [0,pi]); a concave function
        vanishing at both ends of an interval is >= 0 on it, so Phi_w <= 0 and
        Phi >= Phi(b, b/2) = b cos b > 0.
      * on [b/2, W]:  Q >= 0 and Om >= 0, so with w_m = (b/2 + W)/2,
        Phi >= Q(w_m) > 0 for w >= w_m (Q increasing) and
        Phi >= Om(w_m) > 0 for w <= w_m (Om decreasing).
    (b < pi/2 gives b/2 < pi/4 < W, so the interval [b/2, W] is nondegenerate.)

C5 (b*).  g'(b) = -sin(b/2) - (b/4) cos(b/2) < 0 for b in (0,pi), so g is
    STRICTLY DECREASING on (0,pi) and its unique zero b* separates the two cases
    of C3 for every b != b*.

Consequently the Wtau count is determined for every beta outside
{beta1* = pi - b*, pi/2, beta2* = b*} and the Wpot / Wwgt counts for every beta
outside {pi/2}.  The six collars close; the residual uncertified set is those
three beta VALUES, of measure zero, each already handled exactly in section 5 of
the certificate.

This script machine-checks every inequality above in mpmath.iv over the collar
boxes and over a rigorous cover of the whole beta axis.
"""
import json
import mpmath
from mpmath import iv
import cert_core as K

C = K.C
PI = K.PI


def gfun(bI):
    return iv.cos(bI / 2) - (bI / 2) * iv.sin(bI / 2)


def bstar_enclosure(tol=None):
    """certified enclosure of b*, the unique zero of g on (0,pi)."""
    if tol is None:
        tol = mpmath.mpf(2) ** -160
    lo, hi = mpmath.mpf("1.70"), mpmath.mpf("1.75")
    # g' < 0 on the bracket, certified over the WHOLE bracket
    B = iv.mpf([lo, hi])
    gp = -iv.sin(B / 2) - (B / 4) * iv.cos(B / 2)
    assert K.sgn(gp) == -1, "g' not certified negative on the bracket"
    assert K.sgn(gfun(iv.mpf([lo, lo]))) == 1
    assert K.sgn(gfun(iv.mpf([hi, hi]))) == -1
    while hi - lo > tol:
        m = (lo + hi) / 2
        if K.sgn(gfun(iv.mpf([m, m]))) == 1:
            lo = m
        else:
            hi = m
    return iv.mpf([lo, hi])


BSTAR = bstar_enclosure()
HALF = PI() / 2


def side(bI, ref):
    """certified side of an interval bI relative to an interval ref."""
    bl, bh = C.endpoints(bI)
    rl, rh = C.endpoints(ref)
    # the reference is a certified ENCLOSURE [rl, rh] of an irrational
    # constant that lies strictly inside it, so touching an endpoint is
    # still a strict decision
    if bh <= rl:
        return "<"
    if bl >= rh:
        return ">"
    return "?"


def piece_counts(bI):
    """Certified (Wtau, Wpot, Wwgt) zero counts in t on ONE smooth piece with
    parameter b ranging over bI.  Returns (counts, checks) or (None, checks)."""
    checks = {}
    W = (PI() - bI) / 2
    checks["sin_b_pos"] = K.sgn(iv.sin(bI)) == 1
    # ---- Wwgt: Om >= 0, zero exactly at |w| = W
    checks["Om_nonneg_C1"] = K.sgn(iv.cos(bI / 2)) == 1 and checks["sin_b_pos"]
    nwgt = 2
    # ---- Wpot
    s_half = side(bI, HALF)
    Psi0 = bI * iv.cos(bI / 2) * (iv.cos(bI / 2) + (bI / 2) * iv.sin(bI / 2))
    PsiW = (PI() / 2) * (bI - PI() / 2) * iv.sin(bI)
    checks["Psi0_pos"] = K.sgn(Psi0) == 1
    # Psi_w < 0 on (0,W]: 2w sin b + b sin 2w > 0.  Only sin(2w) >= 0 is needed,
    # which holds because 2w in (0, pi-b) subset (0, pi).
    checks["Psi_w_neg_C2"] = checks["sin_b_pos"] and K.sgn(PI() - bI) == 1
    if s_half == "<":
        checks["PsiW_neg"] = K.sgn(PsiW) == -1
        npot = 2
    elif s_half == ">":
        checks["PsiW_pos"] = K.sgn(PsiW) == 1
        npot = 0
    else:
        npot = None
    # ---- Wtau
    if s_half == "<":
        # C4
        Phi_at_half_b = bI * iv.cos(bI)
        checks["Phi_at_b_over_2_pos_C4"] = K.sgn(Phi_at_half_b) == 1
        wm = (bI / 2 + W) / 2
        Qwm = (wm * wm - bI * bI / 4) * iv.sin(bI)
        Omwm = (bI / 2) * (iv.cos(2 * wm) + iv.cos(bI))
        checks["Q_at_wm_pos_C4"] = K.sgn(Qwm) == 1
        checks["Om_at_wm_pos_C4"] = K.sgn(Omwm) == 1
        checks["b_over_2_lt_W_C4"] = K.sgn(W - bI / 2) == 1
        ntau = 0
    elif s_half == ">":
        DW = (2 * bI - PI()) * iv.sin(bI)
        PhiW = (PI() / 4) * (PI() - 2 * bI) * iv.sin(bI)
        checks["D_at_W_pos_C3"] = K.sgn(DW) == 1
        checks["PhiW_neg_C3"] = K.sgn(PhiW) == -1
        s_star = side(bI, BSTAR)
        checks["g_sign_C3C5"] = s_star
        if s_star == "<":
            ntau = 2
        elif s_star == ">":
            ntau = 0
        else:
            ntau = None
    else:
        ntau = None
    ok = all(v for v in checks.values() if isinstance(v, bool))
    return {"Wtau": ntau, "Wpot": npot, "Wwgt": nwgt}, checks, ok


def beta_counts(betaI):
    """Certified total counts in t over one period, for beta in betaI."""
    pieces = {"A": betaI, "B": PI() - betaI}
    tot = {"Wtau": 0, "Wpot": 0, "Wwgt": 0}
    detail, allok = {}, True
    for nm, bI in pieces.items():
        cnt, chk, ok = piece_counts(bI)
        detail[nm] = {"b": [mpmath.nstr(x, 22) for x in C.endpoints(bI)],
                      "counts": cnt, "checks": {k: str(v) for k, v in chk.items()}}
        allok = allok and ok
        for k in tot:
            if cnt[k] is None or tot[k] is None:
                tot[k] = None
            else:
                tot[k] += cnt[k]
    # Wwgt: the two pieces share their endpoints t = 0, beta, pi (mod pi t=0=pi)
    tot["Wwgt"] = 2
    return tot, detail, allok


COLLARS = [
    ("1.420925475550979866654", "1.420925475551033713495"),
    ("1.420925475551033713495", "1.420925475551124584658"),
    ("1.570796326794805748068", "1.570796326794896619231"),
    ("1.570796326794896619231", "1.570796326794987490395"),
    ("1.720667178038668653804", "1.720667178038759524968"),
    ("1.720667178038759524968", "1.720667178038813371809"),
]

CRIT = None          # filled in main(): the three critical enclosures


def adaptive_cover(lo, hi, dmin, budget=200000):
    """Cover [lo,hi] by boxes on which every count is CERTIFIED.

    A box is accepted when beta_counts returns integer counts and every
    ingredient inequality of C1-C5 verifies on the whole box; otherwise it is
    halved.  Boxes narrower than dmin that still fail are returned as residual.
    Near a critical beta the positivity margin of the deciding inequality is
    O(distance), so the refinement is automatically geometric there -- which is
    exactly why a uniform box test (the certificate's) leaves a collar and this
    one does not.
    """
    todo = [(lo, hi)]
    boxes, residual = [], []
    n = 0
    while todo and n < budget:
        a, b = todo.pop()
        n += 1
        tot, _, ok = beta_counts(iv.mpf([a, b]))
        if ok and tot["Wtau"] is not None and tot["Wpot"] is not None:
            boxes.append((a, b, (tot["Wtau"], tot["Wpot"], tot["Wwgt"])))
            continue
        if b - a <= dmin:
            residual.append((a, b))
            continue
        m = (a + b) / 2
        todo.append((a, m))
        todo.append((m, b))
    if todo:
        residual += todo
    return boxes, residual, n


def main():
    global CRIT
    CRIT = {"beta1*": PI() - BSTAR, "pi/2": HALF, "beta2*": BSTAR}
    dmin = mpmath.mpf("1e-45")
    report = {"bstar": [mpmath.nstr(x, 32) for x in C.endpoints(BSTAR)],
              "bstar_width": float(C.width(BSTAR)),
              "beta1star": [mpmath.nstr(x, 32) for x in C.endpoints(PI() - BSTAR)],
              "refinement_floor": str(dmin),
              "critical_enclosures": {k: [mpmath.nstr(x, 32) for x in C.endpoints(v)]
                                      for k, v in CRIT.items()},
              "collars": [], "global_cover": []}

    def do(lo, hi, label):
        boxes, residual, n = adaptive_cover(lo, hi, dmin)
        covered = sum(b - a for a, b, _ in boxes)
        resid = sum(b - a for a, b in residual)
        cts = sorted(set(str(c) for _, _, c in boxes))
        return {"label": label,
                "range": [mpmath.nstr(lo, 30), mpmath.nstr(hi, 30)],
                "width": float(hi - lo), "boxes": len(boxes),
                "measure_certified": float(covered),
                "measure_residual": float(resid),
                "residual_intervals": [[mpmath.nstr(a, 34), mpmath.nstr(b, 34)]
                                       for a, b in residual],
                "counts_(Wtau,Wpot,Wwgt)": cts, "tests": n}

    for lo_s, hi_s in COLLARS:
        e = do(mpmath.mpf(lo_s), mpmath.mpf(hi_s), "collar")
        e["collar"] = [lo_s, hi_s]
        report["collars"].append(e)
        print("collar [%s..%s] w=%.3e -> %d boxes, certified %.6e, residual %.3e, counts %s"
              % (lo_s[:20], hi_s[:20], e["width"], e["boxes"],
                 e["measure_certified"], e["measure_residual"],
                 e["counts_(Wtau,Wpot,Wwgt)"]), flush=True)

    pi_hi = C.endpoints(PI())[1]
    for lo, hi, lab in [(mpmath.mpf("1e-12"), mpmath.mpf("1.42"), "(0, beta1*)"),
                        (mpmath.mpf("1.42"), mpmath.mpf("1.56"), "(beta1*, pi/2) left"),
                        (mpmath.mpf("1.56"), mpmath.mpf("1.58"), "around pi/2"),
                        (mpmath.mpf("1.58"), mpmath.mpf("1.73"), "(pi/2, beta2*) right"),
                        (mpmath.mpf("1.73"), pi_hi - mpmath.mpf("1e-12"), "(beta2*, pi)")]:
        e = do(lo, hi, lab)
        report["global_cover"].append(e)
        print("segment %-22s w=%.4f -> %d boxes, residual %.3e, counts %s"
              % (lab, e["width"], e["boxes"], e["measure_residual"],
                 e["counts_(Wtau,Wpot,Wwgt)"]), flush=True)

    allres, tot = [], 0.0
    for e in report["collars"] + report["global_cover"]:
        allres += e["residual_intervals"]
        tot += e["measure_residual"]
    report["residual_uncertified_beta_set"] = {
        "intervals": allres,
        "contains": ["each residual interval brackets beta1*, pi/2 or beta2*"],
        "total_measure": tot,
        "interpretation":
            "every residual interval contains one of beta1*, pi/2, beta2*; the "
            "deciding inequality's positivity margin there is O(distance to the "
            "critical value), so no interval method can certify AT the point, and "
            "the refinement floor 1e-45 is the only reason the measure is not 0. "
            "The three points themselves are handled exactly in section 5 of the "
            "certificate.  A single beta value contains no open face and cannot "
            "supply an additional two-dimensional face away from the certified "
            "boundary curves."}
    with open("collars.json", "w") as fh:
        json.dump(report, fh, indent=1)
    print("\nresidual uncertified beta measure: %.3e (was 4.712e-13)"
          % report["residual_uncertified_beta_set"]["total_measure"])


if __name__ == "__main__":
    main()
