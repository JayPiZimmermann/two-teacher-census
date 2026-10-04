"""Independent falsification gate for the normalized selected-wall forms.

The DFS evaluates algebraically normalized Taylor forms, while the reference
side here calls the literal angle rows and the independently transcribed
``kernel_schur`` selector.  At thin points the exact forced powers are

    locus = D^6 * locus_bar,
    kernel_T00 = D^10 * T00_bar,
    kernel_det = D^20 * det_bar.

The gate checks both thin points and deterministic interior probes of interval
boxes.  It is a falsifier, not the certificate; replay still recomputes every
leaf with directed interval arithmetic.
"""
import math
import random

import census_cert as X
import minimum_wall_dfs as W
import noncentered as J


def contains(interval, value):
    return J.lo(interval) <= J.lo(value) and J.hi(value) <= J.hi(interval)


def overlaps(left, right):
    return J.lo(left) <= J.hi(right) and J.lo(right) <= J.hi(left)


def reference(beta, s, gap):
    forms = W.wall_forms_literal(X.thin(beta), X.thin(s), X.thin(gap))
    d = X.thin(gap)
    return forms[0] / d ** 6, forms[1] / d ** 10, forms[2] / d ** 20


def smooth(beta, s, gap, pad=0.0):
    angles = (s, s + gap, s - beta, s + gap - beta)
    return min(abs(math.sin(z)) for z in angles) > pad


def main():
    rng = random.Random(20260810)
    thin_fail = 0
    box_fail = 0
    thin_seen = 0
    box_seen = 0
    worst_relative_width = [0.0, 0.0, 0.0]

    while thin_seen < 500:
        beta = rng.uniform(0.08, 3.05)
        s = rng.uniform(0.05, 2 * math.pi - 0.05)
        gap = 10 ** rng.uniform(-4.0, math.log10(0.7))
        if not smooth(beta, s, gap, 2e-5):
            continue
        got = W.wall_forms_normalized(
            X.thin(beta), X.thin(s), X.thin(gap))
        if got is None:
            continue
        want = reference(beta, s, gap)
        thin_seen += 1
        for i in range(3):
            # At tiny gaps the independently evaluated literal expression,
            # followed by division by D^20, has a wider rounding enclosure
            # than the Taylor form.  Both are enclosures of the same scalar;
            # requiring overlap is the meaningful independent thin-point
            # falsifier.  Interior point references below must be contained
            # by the genuinely non-thin box enclosure.
            if not overlaps(got[i], want[i]):
                thin_fail += 1
            scale = max(1.0, abs(float(J.mid(want[i]))))
            width = float(J.hi(got[i]) - J.lo(got[i])) / scale
            worst_relative_width[i] = max(worst_relative_width[i], width)

    while box_seen < 120:
        beta = rng.uniform(0.12, 3.0)
        s = rng.uniform(0.12, 2 * math.pi - 0.12)
        gap = 10 ** rng.uniform(-3.5, math.log10(0.65))
        rb = min(2e-5, beta / 20)
        rs = min(2e-5, gap / 20)
        rd = min(1e-5, gap / 10)
        if not smooth(beta, s, gap, 8 * max(rb, rs, rd)):
            continue
        B = W.M.mk(beta - rb, beta + rb)
        S = W.M.mk(s - rs, s + rs)
        D = W.M.mk(gap - rd, gap + rd)
        got = W.wall_forms_normalized(B, S, D)
        if got is None:
            continue
        box_seen += 1
        for k in range(9):
            qbeta = beta + rb * (2 * ((3 * k + 1) % 11) / 10 - 1)
            qs = s + rs * (2 * ((5 * k + 2) % 11) / 10 - 1)
            qgap = gap + rd * (2 * ((7 * k + 3) % 11) / 10 - 1)
            want = reference(qbeta, qs, qgap)
            for i in range(3):
                if not contains(got[i], want[i]):
                    box_fail += 1

    print("normalized wall gate: thin=%d failures=%d boxes=%d "
          "interior_failures=%d" % (thin_seen, thin_fail, box_seen, box_fail))
    print("worst thin relative widths: locus %.3e T00 %.3e det %.3e" %
          tuple(worst_relative_width))
    if thin_fail or box_fail:
        print("GATE FAILED")
        return 1
    print("GATE PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
