"""PROSPECTING (not evidence): WHICH STAGE loses the small-gap SCHUR test?

The published reading of the small-gap wall is that `|det|/|det'|` -- "the
distance to the nearest zero of the determinant" -- falls to `1e-25` at
`D = 0.032`, i.e. that the quantity itself is degenerate and no enclosure form
can help.  An independent central-difference measurement (`row7a_band.py`)
puts the linearised distance to a zero at `1.0e-3` there, not `1e-25`: `|det|`
falls 31 orders as the gap shrinks and so does its DERIVATIVE, because the
whole determinant is a degree-four form in a kernel vector that is itself going
to zero.  The two readings cannot both be right, so this script prices the
enclosure stage by stage.

For one family it reports, at a thin point and at boxes of decreasing side, the
RELATIVE enclosure width of every intermediate quantity of
`kernel_schur.kernel_teacher_schur_detT`:

    S0, S1  (the kernel vector -- `num` at two arguments)
    N0, N1  (the Cramer entries)
    T00, T11, T01, det

A relative width that is already large at the THIN point is a cancellation in
the evaluation and no box shrinking can fix it.  A relative width that falls
linearly with the box side is dependency, and dependency is what a mean-value
or a rewritten form attacks.  The two have completely different repairs, and
the published prose names neither.

Run: `python3 row7a_stagewidth.py [bits]`
"""
import sys

import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS
import kernel_schur_mv as MV


def relw(x):
    m = abs(J.mid(x))
    if m == 0:
        return mpmath.inf
    return J.width(x) / m


def stages(B, T0, T1):
    """Every intermediate of the determinant, in evaluation order."""
    D = T0 - T1
    AD = KS._atoms(D)
    phiD, hD = AD.phi, AD.h
    M = X.PIv * X.PIv - phiD * phiD
    S0 = -KS.num(D, T0 - B)
    S1 = KS.num(D, T0)

    def V(g):
        return S0 * KS._atoms(g).phi + S1 * KS._atoms(g - B).phi

    def L(g):
        return S0 * KS._abs_sin(g) + S1 * KS._abs_sin(g - B)

    V0, V1 = V(T0), V(T1)
    N0 = X.PIv * V0 - phiD * V1
    N1 = X.PIv * V1 - phiD * V0
    K = 2 * KS._abs_sin(D) - phiD
    base = M * N0 * N1 * K
    T00 = base - M * M * N0 * (2 * L(T0) - V0) - X.PIv * hD * hD * N0 * N0
    T11 = base - M * M * N1 * (2 * L(T1) - V1) - X.PIv * hD * hD * N1 * N1
    T01 = -base - phiD * hD * hD * N0 * N1
    det = T00 * T11 - T01 * T01
    return [("M", M), ("S0", S0), ("S1", S1), ("V0", V0), ("V1", V1),
            ("N0", N0), ("N1", N1), ("T00", T00), ("T11", T11),
            ("T01", T01), ("det", det)]


FAMS = [
    ("F3 (3.08,-0.08)", 3.08, 6.220945, 0.0315),
    ("F7 (3.08,-0.92)", 3.08, 6.220945, 0.0315),
    ("F2 (0.75,-0.02)", 0.75, 5.671000, 0.2019),
]


def detval(b, s, d):
    """Thin-point value of the determinant (160-bit interval midpoint)."""
    v, _ = KS.kernel_teacher_schur_detT(
        X.iv.mpf([b, b]), X.iv.mpf([s + d, s + d]), X.iv.mpf([s, s]))
    return J.mid(v)


def true_range(b, s, d, h, n=9):
    """The TRUE range of the determinant over the box, by dense sampling.

    Sampling under-estimates a range (conventions §1p), so this is a LOWER
    bound on the true range and therefore an UPPER bound on how much an ideal
    enclosure could gain -- which is the direction that keeps the headroom
    claim honest."""
    lo = hi = None
    for i in range(n):
        for k in range(n):
            for j in range(3):
                v = detval(b - h / 2 + h * j / 2.0,
                           s - h / 2 + h * i / (n - 1.0),
                           d - h / 2 + h * k / (n - 1.0))
                lo = v if lo is None else min(lo, v)
                hi = v if hi is None else max(hi, v)
    return lo, hi


def headroom(allfam=False, side="1e-4", n=5):
    """How much of the enclosure width is DEPENDENCY, at the families where
    the sign test fails: enclosure width against the sampled true range.

    With `allfam` the same ratio is reported for EVERY certified family at one
    box side, which is what answers "is the enclosure's growth constant a
    measurable property of the family?"  It is: it tracks the student gap."""
    import kernel_schur_headroom as H
    fam = [f for f in H.families() if not f[6]]
    fam.sort(key=lambda f: float(f[4]) - float(f[3]))
    if allfam:
        h = mpmath.mpf(side)
        print("\nDEPENDENCY FACTOR at box side %s, all families, by gap"
              % side)
        print("  %-24s %8s %12s %12s %12s %10s"
              % ("witness", "D", "|det|", "true range", "mv width", "factor"))
        for (tag, beta, y, th1, th0, _s, _f) in fam:
            b = mpmath.mpf(repr(beta))
            ss = mpmath.mpf(repr(th1))
            d = mpmath.mpf(repr(th0)) - ss
            lo, hi = true_range(b, ss, d, h, n=n)
            rng = hi - lo
            B = X.iv.mpf([b - h / 2, b + h / 2])
            mv = MV.detT_meanvalue(B, ss - h / 2, ss + h / 2,
                                   d - h / 2, d + h / 2)
            w = J.width(mv)
            v = detval(b, ss, d)
            fac = w / rng if rng > 0 else mpmath.inf
            print("  %-24s %8.4f %12.3e %12.3e %12.3e %10.2e"
                  % (tag.replace("replay_C_", "").replace(".log", ""),
                     float(d), float(abs(v)), float(rng), float(w),
                     float(fac)))
        return
    print("\nDEPENDENCY HEADROOM: enclosure width / sampled true range")
    print("  %-24s %8s %6s %12s %12s %12s %10s"
          % ("witness", "D", "side", "|det|", "true range", "mv width",
             "factor"))
    for (tag, beta, y, th1, th0, _s, _f) in fam[:3] + [fam[len(fam) // 2]]:
        b = mpmath.mpf(repr(beta))
        s = mpmath.mpf(repr(th1))
        d = mpmath.mpf(repr(th0)) - s
        for h in (mpmath.mpf("1e-4"), mpmath.mpf("1e-6")):
            lo, hi = true_range(b, s, d, h)
            rng = hi - lo
            B = X.iv.mpf([b - h / 2, b + h / 2])
            mv = MV.detT_meanvalue(B, s - h / 2, s + h / 2,
                                   d - h / 2, d + h / 2)
            w = J.width(mv)
            v = detval(b, s, d)
            fac = w / rng if rng > 0 else mpmath.inf
            print("  %-24s %8.4f %6.0e %12.3e %12.3e %12.3e %10.2e"
                  % (tag.replace("replay_C_", "").replace(".log", ""),
                     float(d), float(h), float(abs(v)), float(rng), float(w),
                     float(fac)))


def main():
    bits = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    J.set_prec(bits)
    if len(sys.argv) > 2 and sys.argv[2] == "all":
        headroom(allfam=True)
        return 0
    print("working precision: %d bits" % bits)
    # the family midpoints are re-derived from the replayed witnesses, so the
    # exact s below matters only as a place to stand: any point of the family
    # curve has the same stage structure.
    import kernel_schur_headroom as H
    fam = [f for f in H.families() if not f[6]]
    fam.sort(key=lambda f: float(f[4]) - float(f[3]))
    picks = fam[:3] + [fam[len(fam) // 2]] + [fam[-1]]
    for (tag, beta, y, th1, th0, _s, _f) in picks:
        b = mpmath.mpf(repr(beta))
        s = mpmath.mpf(repr(th1))
        d = mpmath.mpf(repr(th0)) - s
        print("\n=== %s beta=%.2f y=%.2f D=%.4f" %
              (tag.replace("replay_C_", "").replace(".log", ""), beta, y,
               float(d)))
        print("  %-5s %12s | %s" % ("side", "|det|",
                                    "relative enclosure width per stage"))
        for h in (mpmath.mpf(0), mpmath.mpf("1e-4"), mpmath.mpf("1e-6"),
                  mpmath.mpf("1e-8"), mpmath.mpf("1e-10")):
            B = X.iv.mpf([b - h / 2, b + h / 2])
            T1 = X.iv.mpf([s - h / 2, s + h / 2])
            Dm = X.iv.mpf([d - h / 2, d + h / 2])
            st = stages(B, T1 + Dm, T1)
            det = st[-1][1]
            txt = "  ".join("%s %.1e" % (n, float(relw(v))) for n, v in st
                            if n in ("S0", "S1", "N0", "T00", "T01", "det"))
            mv = MV.detT_meanvalue(B, s - h / 2, s + h / 2,
                                   d - h / 2, d + h / 2)
            print("  %-5s %12.3e | %s | mv relw %.1e  fires lit=%s mv=%s"
                  % ("thin" if h == 0 else ("%.0e" % float(h)),
                     float(abs(J.mid(det))), txt, float(relw(mv)),
                     X.sgn(det) != 0, X.sgn(mv) != 0))
    headroom()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
