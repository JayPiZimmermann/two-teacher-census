"""GATE for the mean-value enclosure of the mass-free Schur determinant.

Three checks, all before the form is allowed to discharge a box:

  (a) VALUE.  The AD pass reproduces the literal determinant of
      `kernel_schur.kernel_teacher_schur_detT` at points.  If the two disagree
      the AD rebuild is not the same polynomial.
  (b) DERIVATIVE.  The AD derivatives `d det/ds` and `d det/dD` agree with
      central differences of the literal determinant.  This certificate has
      already paid once for a derivative enclosure that was merely invalid, so
      it is checked rather than trusted.
  (c) ENCLOSURE.  On boxes, the mean-value form must CONTAIN the literal
      enclosure's true value -- checked by sampling the literal determinant
      inside the box and confirming every sample lies in the mean-value
      interval.  A tighter enclosure that fails to contain is worthless.
  (d) ABSOLUTE-SINE TRUTH.  The literal evaluator agrees with an independent
      rebuild that obtains `|sin|` from `noncentered.Atoms.asin`.  This catches
      a missing odd-half-branch orientation in both a literal evaluator and a
      self-consistent AD copy.

Run: `python3 kernel_schur_mv_gate.py`
"""
import random

import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS
import kernel_schur_mv as MV

random.seed(20260808)


def pt(v):
    m = v if isinstance(v, mpmath.mpf) else mpmath.mpf(repr(v))
    return X.iv.mpf([m, m])


def truth(beta, T0, T1):
    """Independent selector rebuild using the canonical atom's `.asin`."""
    Dm = T0 - T1
    AD, _ = J.atoms(Dm)
    mass_det = X.PIv * X.PIv - AD.phi * AD.phi
    s0 = -KS.num(Dm, T0 - beta)
    s1 = KS.num(Dm, T0)

    def load(T):
        A0, _ = J.atoms(T)
        A1, _ = J.atoms(T - beta)
        return (s0 * A0.phi + s1 * A1.phi,
                s0 * A0.asin + s1 * A1.asin)

    P0, L0 = load(T0)
    P1, L1 = load(T1)
    n0 = X.PIv * P0 - AD.phi * P1
    n1 = X.PIv * P1 - AD.phi * P0
    band = 2 * AD.asin - AD.phi
    base = mass_det * n0 * n1 * band
    t00 = base - mass_det * mass_det * n0 * (2 * L0 - P0) \
        - X.PIv * AD.h * AD.h * n0 * n0
    t11 = base - mass_det * mass_det * n1 * (2 * L1 - P1) \
        - X.PIv * AD.h * AD.h * n1 * n1
    t01 = -base - AD.phi * AD.h * AD.h * n0 * n1
    return t00 * t11 - t01 * t01, (t00, t11, t01)


def main():
    wv = wd = wt = mpmath.mpf(0)
    n = 0
    for _ in range(200):
        b = random.uniform(0.1, 3.05)
        s = random.uniform(0.2, 6.0)
        d = random.uniform(0.05, 3.10)
        B, T1, Dm = pt(b), pt(s), pt(d)
        T0 = T1 + Dm
        lit, parts = KS.kernel_teacher_schur_detT(B, T0, T1)
        direct, direct_parts = truth(B, T0, T1)
        for a, c in zip((lit,) + parts[:3], (direct,) + direct_parts):
            scale = max(abs(J.mid(c)), mpmath.mpf("1e-30"))
            wt = max(wt, abs(J.mid(a) - J.mid(c)) / scale)
        v, gs, gd = MV.detT_D2(B, T0, T1, Dm)
        sc = max(abs(J.mid(lit)), mpmath.mpf("1e-30"))
        wv = max(wv, abs(J.mid(v) - J.mid(lit)) / sc)
        e = mpmath.mpf("1e-18")
        for coord, got in ((0, gs), (1, gd)):
            if coord == 0:
                a, _ = KS.kernel_teacher_schur_detT(B, pt(s + e) + Dm, pt(s + e))
                c, _ = KS.kernel_teacher_schur_detT(B, pt(s - e) + Dm, pt(s - e))
            else:
                a, _ = KS.kernel_teacher_schur_detT(B, T1 + pt(d + e), T1)
                c, _ = KS.kernel_teacher_schur_detT(B, T1 + pt(d - e), T1)
            fd = (J.mid(a) - J.mid(c)) / (2 * e)
            sc2 = max(abs(fd), mpmath.mpf("1e-18"))
            wd = max(wd, abs(J.mid(got) - fd) / sc2)
        n += 1
    print("(a) AD value vs literal determinant : %d samples, worst relative "
          "residual %.3e" % (n, float(wv)))
    print("(b) AD derivative vs central differences : worst relative residual "
          "%.3e" % float(wd))
    if wv > mpmath.mpf("1e-25") or wd > mpmath.mpf("1e-8"):
        print("GATE FAILED")
        return 1

    bad = 0
    tested = 0
    for _ in range(120):
        b = random.uniform(0.3, 2.8)
        s = random.uniform(0.3, 5.8)
        d = random.uniform(0.1, 3.0)
        w = random.choice([1e-2, 1e-3, 1e-4])
        B = pt(b)
        mv = MV.detT_meanvalue(B, mpmath.mpf(repr(s)), mpmath.mpf(repr(s + w)),
                               mpmath.mpf(repr(d)), mpmath.mpf(repr(d + w)))
        lo, hi = J.endpoints(mv)
        for i in range(5):
            for k in range(5):
                ss = pt(s + w * i / 4.0)
                dd = pt(d + w * k / 4.0)
                lit, _ = KS.kernel_teacher_schur_detT(B, ss + dd, ss)
                m = J.mid(lit)
                tested += 1
                if m < lo or m > hi:
                    bad += 1
    print("(c) mean-value form contains the true value : %d samples, %d "
          "outside" % (tested, bad))
    print("(d) literal vs canonical `.asin` rebuild : worst relative residual "
          "%.3e" % float(wt))
    if bad or wt > mpmath.mpf("1e-25"):
        print("GATE FAILED: the mean-value form is not an enclosure")
        return 1
    print("all gates passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
