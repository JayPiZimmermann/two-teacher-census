"""GATE for `census_cert.sep_res_rowform` (the division-free row form).

Two things are checked, both before the row form is allowed to exclude a box:

1. IDENTITY.  `E_i = massDet(D) * G_i` exactly, where `G_i` is what `sep_res`
   returns.  Checked at POINT boxes so both sides are essentially exact, at
   random teachers and student pairs, in both mass sectors.  A failure here
   means the row form is not the same equation and must not be used.

2. TIGHTNESS.  On genuine (non-degenerate) boxes, the ratio of the row form's
   enclosure width to the width of `massDet * G_i` computed the Cramer way.
   The row form never divides by `massDet`, so it should be at least as tight;
   this reports by how much, which is the whole reason for the second cut.

Nothing here certifies anything: it is a falsifier for the identity and a
measurement of the width, and it is prospecting in the artifact rule's sense.
Run: `python3 validate_rowform.py`
"""
import random

import mpmath

import census_cert as X
import noncentered as J

random.seed(20260808)


def point(v):
    m = v if isinstance(v, mpmath.mpf) else mpmath.mpf(repr(v))
    return X.iv.mpf([m, m])


def box(lo, hi):
    return X.iv.mpf([mpmath.mpf(repr(lo)), mpmath.mpf(repr(hi))])


def jac_gates():
    """The DERIVATIVE gates for `sep_res_rowform_jac`, two of them.

    (a) against the audited Cramer Jacobian, through the exact relations
        dE/ds = M * dG/ds  and  dE_i/dD = (dM/dD) * G_i + M * dG_i/dD;
    (b) against central finite differences at the same point.

    An invalid derivative enclosure is the failure that voided verdicts in
    this certificate once already, so it is gated before use, not after.
    """
    worst_rel = mpmath.mpf(0)
    worst_fd = mpmath.mpf(0)
    n = 0
    for _ in range(300):
        beta = random.uniform(0.05, 3.09)
        y = random.uniform(-0.99, 0.99)
        if abs(y) < 0.02:
            continue
        th1 = random.uniform(0.0, 6.28)
        D = random.uniform(0.05, 3.05)
        s0, s1 = X.masses_at(point(y))
        Tc = X.Teacher(point(beta), s0, s1)
        T1, Dm = point(th1), point(D)
        T0 = T1 + Dm
        r = X.sep_res(Tc, T0, T1, None, Dm=Dm, want_jac=True)
        if r is None:
            continue
        G0, G1, M = r[0], r[1], r[4]
        hD, phiD = r[5], r[6]
        j00, j01, j10, j11 = r[8]
        E0, E1, ds0, dd0, ds1, dd1 = X.sep_res_rowform_jac(
            Tc, T0, T1, None, Dm=Dm)
        q = 2 * phiD * hD
        for got, want in ((ds0, M * (j00 + j01)),
                          (ds1, M * (j10 + j11)),
                          (dd0, q * G0 + M * j00),
                          (dd1, q * G1 + M * j10)):
            d = J.mid(got) - J.mid(want)
            sc = max(abs(J.mid(want)), mpmath.mpf("1e-30"))
            worst_rel = max(worst_rel, abs(d) / sc)
        # central differences in each coordinate
        e = mpmath.mpf("1e-20")
        for coord in (0, 1):
            if coord == 0:
                a = X.sep_res_rowform(Tc, point(th1 + e) + Dm,
                                      point(th1 + e), None, Dm=Dm)
                b = X.sep_res_rowform(Tc, point(th1 - e) + Dm,
                                      point(th1 - e), None, Dm=Dm)
                got = (ds0, ds1)
            else:
                a = X.sep_res_rowform(Tc, T1 + point(D + e), T1, None,
                                      Dm=point(D + e))
                b = X.sep_res_rowform(Tc, T1 + point(D - e), T1, None,
                                      Dm=point(D - e))
                got = (dd0, dd1)
            for i in (0, 1):
                fd = (J.mid(a[i]) - J.mid(b[i])) / (2 * e)
                sc = max(abs(fd), mpmath.mpf("1e-20"))
                worst_fd = max(worst_fd, abs(J.mid(got[i]) - fd) / sc)
        n += 1
    print("derivative vs audited Cramer Jacobian : %d samples, worst relative "
          "residual %.3e" % (n, float(worst_rel)))
    print("derivative vs central differences     : worst relative residual "
          "%.3e" % float(worst_fd))
    return 0 if (worst_rel < mpmath.mpf("1e-20")
                 and worst_fd < mpmath.mpf("1e-10")) else 1


def main():
    # ---- 1. the identity, at point boxes
    worst = mpmath.mpf(0)
    n_ok = 0
    for _ in range(400):
        beta = random.uniform(0.05, 3.09)
        y = random.uniform(-0.99, 0.99)
        if abs(y) < 0.02:
            continue
        th1 = random.uniform(0.0, 6.28)
        D = random.uniform(0.02, 3.10)
        s0, s1 = X.masses_at(point(y))
        Tc = X.Teacher(point(beta), s0, s1)
        T1, Dm = point(th1), point(D)
        T0 = T1 + Dm
        r = X.sep_res(Tc, T0, T1, None, Dm=Dm)
        if r is None:
            continue
        G0, G1, _, _, M = r[0], r[1], r[2], r[3], r[4]
        E0, E1 = X.sep_res_rowform(Tc, T0, T1, None, Dm=Dm)
        for E, G in ((E0, G0), (E1, G1)):
            d = J.mid(E) - J.mid(M) * J.mid(G)
            scale = max(abs(J.mid(E)), mpmath.mpf("1e-30"))
            worst = max(worst, abs(d) / scale)
        n_ok += 1
    print("identity  E_i = massDet * G_i : %d point samples, worst relative "
          "residual %.3e" % (n_ok, float(worst)))
    if worst > mpmath.mpf("1e-20"):
        print("GATE FAILED: the row form is NOT the same equation")
        return 1

    # ---- 2. the width, on genuine boxes
    ratios = []
    nodiv = 0
    for _ in range(400):
        beta = random.uniform(0.05, 3.09)
        y = random.uniform(-0.99, 0.99)
        if abs(y) < 0.02:
            continue
        th1 = random.uniform(0.0, 6.0)
        D = random.uniform(0.05, 3.05)
        w = random.choice([1e-1, 1e-2, 1e-3, 1e-4])
        s0, s1 = X.masses_at(point(y))
        Tc = X.Teacher(point(beta), s0, s1)
        T1 = box(th1, th1 + w)
        Dm = box(D, D + w)
        T0 = T1 + Dm
        E0, E1 = X.sep_res_rowform(Tc, T0, T1, None, Dm=Dm)
        r = X.sep_res(Tc, T0, T1, None, Dm=Dm)
        if r is None:
            nodiv += 1
            continue
        G0, G1, M = r[0], r[1], r[4]
        for E, G in ((E0, G0), (E1, G1)):
            wr, wc = J.width(E), J.width(M * G)
            if wc > 0:
                ratios.append(float(wr / wc))
    ratios.sort()
    n = len(ratios)
    print("width ratio  row / (massDet * Cramer) : n = %d, median %.3f, "
          "p90 %.3f, max %.3f" % (n, ratios[n // 2], ratios[int(0.9 * n)],
                                  ratios[-1]))
    print("boxes where the Cramer form gives up (massDet straddles 0) but the "
          "row form still evaluates: %d" % nodiv)
    return jac_gates()


if __name__ == "__main__":
    raise SystemExit(main())
