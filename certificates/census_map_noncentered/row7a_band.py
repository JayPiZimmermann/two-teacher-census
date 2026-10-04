"""PROSPECTING (not evidence): where can the row-7a SCHUR test fire at all?

Re-measures, independently of the earlier runs, the two numbers that decide
the row-7a exclusion sweep:

  (1) the REQUIRED BOX SIDE per certified family -- the largest box side `h`
      at which the enclosure of `generalJKernelTeacherSchurDetT` over the box
      still misses zero -- for the LITERAL enclosure and for the MEAN-VALUE
      form, and, unlike the earlier table, ALSO with `beta` carried as an
      interval of width `h` rather than as a thin value.  A sweep over a face
      subdivides `beta`; a required side measured at thin `beta` is a
      statement about a different sweep.

  (2) the linearised DISTANCE TO A ZERO of the determinant,
      `|det| / (|d/ds| + |d/dD| + |d/dbeta|)`, with the derivatives taken by
      CENTRAL DIFFERENCES of the literal determinant at 160-bit precision --
      deliberately NOT through the forward-mode AD path of
      `kernel_schur_mv.py`, so that the two agree only if the mathematics
      agrees and not because they share code (conventions skill 0f).

VALIDITY DOMAIN.  Point/box evaluations at the family midpoints recorded in
the replayed witnesses; 160-bit interval arithmetic; central differences at
step 1e-20 with mpmath.mp at 180 bits.  Nothing here is a certificate: it
prices a sweep and characterises a boundary.

Run: `python3 row7a_band.py`
"""
import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS
import kernel_schur_mv as MV
import kernel_schur_headroom as H

SIDES = [mpmath.mpf("1e-1"), mpmath.mpf("1e-2"), mpmath.mpf("1e-3"),
         mpmath.mpf("1e-4"), mpmath.mpf("1e-5"), mpmath.mpf("1e-6"),
         mpmath.mpf("1e-7"), mpmath.mpf("1e-8")]


def pt(v):
    m = v if isinstance(v, mpmath.mpf) else mpmath.mpf(repr(v))
    return X.iv.mpf([m, m])


def box(c, h):
    return X.iv.mpf([c - h / 2, c + h / 2])


def fires_literal(B, s, d, h, wide_beta):
    Bb = box(J.mid(B), h) if wide_beta else B
    T1 = box(s, h)
    Dm = box(d, h)
    try:
        v, _ = KS.kernel_teacher_schur_detT(Bb, T1 + Dm, T1)
    except Exception:
        return False
    return X.sgn(v) != 0


def fires_mv(B, s, d, h, wide_beta):
    Bb = box(J.mid(B), h) if wide_beta else B
    try:
        v = MV.detT_meanvalue(Bb, s - h / 2, s + h / 2, d - h / 2, d + h / 2)
    except Exception:
        return False
    return X.sgn(v) != 0


def required(fn, B, s, d, wide_beta):
    """Largest side in SIDES at which the test fires; None if none does."""
    for h in SIDES:
        if fn(B, s, d, h, wide_beta):
            return h
    return None


def detval(b, s, d):
    v, _ = KS.kernel_teacher_schur_detT(pt(b), pt(s + d), pt(s))
    return J.mid(v)


def main():
    fam = [f for f in H.families() if not f[6]]
    print("non-exact-fit certified families read from replayed witnesses: %d"
          % len(fam))
    e = mpmath.mpf("1e-20")
    rows = []
    for (tag, beta, y, th1, th0, schur, _fit) in fam:
        b = mpmath.mpf(repr(beta))
        s = mpmath.mpf(repr(th1))
        d = mpmath.mpf(repr(th0)) - s
        B = pt(b)
        v = detval(b, s, d)
        gs = (detval(b, s + e, d) - detval(b, s - e, d)) / (2 * e)
        gd = (detval(b, s, d + e) - detval(b, s, d - e)) / (2 * e)
        gb = (detval(b + e, s, d) - detval(b - e, s, d)) / (2 * e)
        gsum = abs(gs) + abs(gd) + abs(gb)
        dist = abs(v) / gsum if gsum > 0 else mpmath.inf
        rows.append({
            "tag": tag, "beta": beta, "y": y, "D": float(d), "s": float(s),
            "det": v, "dist": dist,
            "lit_thin": required(fires_literal, B, s, d, False),
            "mv_thin": required(fires_mv, B, s, d, False),
            "lit_wide": required(fires_literal, B, s, d, True),
            "mv_wide": required(fires_mv, B, s, d, True),
        })

    def tab(key):
        buckets = {}
        for r in rows:
            k = "unresolved" if r[key] is None else ("%.0e" % float(r[key]))
            buckets[k] = buckets.get(k, 0) + 1
        return buckets

    print("\nREQUIRED BOX SIDE, counts by decade")
    print("  %-28s %s" % ("form", "buckets"))
    for key, name in (("lit_thin", "literal, beta THIN"),
                      ("mv_thin", "mean-value, beta THIN"),
                      ("lit_wide", "literal, beta WIDE (width h)"),
                      ("mv_wide", "mean-value, beta WIDE (width h)")):
        b = tab(key)
        order = sorted([k for k in b if k != "unresolved"],
                       key=lambda z: -float(z))
        txt = "  ".join("%s:%d" % (k, b[k]) for k in order)
        if "unresolved" in b:
            txt += "   UNRESOLVED(<1e-8):%d" % b["unresolved"]
        print("  %-28s %s" % (name, txt))

    rows.sort(key=lambda r: r["D"])
    print("\nPER FAMILY, sorted by student gap D "
          "(dist = |det| / (|d/ds|+|d/dD|+|d/dbeta|), central differences)")
    print("  %-24s %6s %6s %8s %11s %10s %8s %8s %8s %8s"
          % ("witness", "beta", "y", "D", "|det|", "dist", "lit/thin",
             "mv/thin", "lit/wide", "mv/wide"))
    for r in rows:
        def f(k):
            return "--" if r[k] is None else ("%.0e" % float(r[k]))
        print("  %-24s %6.2f %6.2f %8.4f %11.3e %10.3e %8s %8s %8s %8s"
              % (r["tag"].replace("replay_C_", "").replace(".log", ""),
                 r["beta"], r["y"], r["D"], float(abs(r["det"])),
                 float(r["dist"]), f("lit_thin"), f("mv_thin"),
                 f("lit_wide"), f("mv_wide")))

    print("\nD-SCALING of the distance-to-a-zero (this is what closes the "
          "small-gap band):")
    for lo, hi in ((0.0, 0.05), (0.05, 0.15), (0.15, 0.3), (0.3, 0.5),
                   (0.5, 1.0), (1.0, 10.0)):
        sel = [r for r in rows if lo <= r["D"] < hi]
        if not sel:
            continue
        ds = [float(r["dist"]) for r in sel]
        print("  D in [%.2f, %.2f): %2d families, dist min %.2e max %.2e"
              % (lo, hi, len(sel), min(ds), max(ds)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
