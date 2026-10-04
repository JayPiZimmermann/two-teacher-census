"""Numerical falsification gate for the Cramer selector transcription.

Checks at deterministic off-kink samples that

* the division-free formulas equal ``massDet^2 * T00`` and
  ``massDet^4 * detT`` computed by the older mass-carrying Schur code;
* all four interval-AD partials agree with centered finite differences; and
* the first-order enclosure contains direct thin evaluations at deterministic
  interior probes of random boxes.
* the absolute and teacher-relative chart transcriptions agree; and
* the small-gap Taylor row enclosure contains the original row divided by
  ``D^3``.

This gate is not a proof and is not part of a certificate.  It is a required
falsifier before emitting certificates with ``minimum_dfs.py``.
"""
import math
import random

import census_cert as X
import minimum_dfs as D
import noncentered as J


def midpoint(z):
    return float(J.mid(z))


def mass_carrying(beta, y, s, gap):
    B, Y, S, G = map(X.thin, (beta, y, s, gap))
    s0, s1 = X.masses_at(Y)
    teacher = X.Teacher(B, s0, s1)
    T0 = S + G
    residual = X.sep_res(teacher, T0, S, None, Dm=G)
    phi, h, slope = None, None, None
    atom_d, _ = J.atoms(G)
    phi, h, slope = atom_d.phi, atom_d.h, atom_d.sA
    P0 = teacher.load(T0)[0]
    P1 = teacher.load(S)[0]
    a0, _ = J.atoms(T0)
    a0b, _ = J.atoms(T0 - B)
    a1, _ = J.atoms(S)
    a1b, _ = J.atoms(S - B)
    L0 = s0 * a0.asin + s1 * a0b.asin
    L1 = s0 * a1.asin + s1 * a1b.asin
    c0, c1, mass_det = residual[2], residual[3], residual[4]
    raw00 = -(c0 * c1 * slope) + c0 * (P0 - 2 * L0)
    raw11 = -(c0 * c1 * slope) + c1 * (P1 - 2 * L1)
    raw01 = c0 * c1 * slope
    t00 = mass_det * raw00 - X.PIv * h * h * c0 * c0
    t11 = mass_det * raw11 - X.PIv * h * h * c1 * c1
    t01 = mass_det * raw01 - phi * h * h * c0 * c1
    det = t00 * t11 - t01 * t01
    return midpoint(mass_det * mass_det * t00), midpoint(mass_det ** 4 * det)


def selectors(point):
    args = tuple(X.thin(v) for v in point)
    out = D.cramer_selectors_d4(*args)
    return [midpoint(out[i].v) for i in (0, 3)], out


def selectors_mode(point, relative):
    args = tuple(X.thin(v) for v in point)
    out = D.cramer_selectors_d4(*args, relative)
    return [midpoint(out[i].v) for i in (0, 3)], out


def rows_mode(point, relative):
    args = tuple(X.thin(v) for v in point)
    out = D.coordinate_rowform_d4(*args, relative)
    return [midpoint(out[i].v) for i in (0, 1)], out


def contains(enclosure, value):
    return J.lo(enclosure) <= value <= J.hi(enclosure)


def main():
    rng = random.Random(20260810)
    formula_err = 0.0
    deriv_err = 0.0
    contain_fail = 0
    literal_err = 0.0
    chart_err = 0.0
    row_deriv_err = 0.0
    row_contain_fail = 0
    kraw_contraction_fail = 0
    collar_contain_fail = 0
    quotient_contain_fail = 0
    samples = []
    while len(samples) < 120:
        p = (rng.uniform(0.12, 3.02), rng.uniform(-0.98, 0.92),
             rng.uniform(0.2, 6.0), rng.uniform(0.08, 3.02))
        # Stay away from atom kinks, where centered finite differences are a
        # poor falsifier for a one-sided derivative enclosure.
        atoms = (p[2], p[2] + p[3], p[2] - p[0], p[2] + p[3] - p[0], p[3])
        if min(abs(math.sin(t)) for t in atoms) < 0.06:
            continue
        samples.append(p)

    for point in samples:
        direct, ad = selectors(point)
        literal = D.cramer_selectors(*(X.thin(v) for v in point))
        for got, enclosed in zip(direct, (literal[0], literal[3])):
            other = midpoint(enclosed)
            literal_err = max(literal_err,
                              abs(got - other) / max(1.0, abs(got), abs(other)))
        old = mass_carrying(*point)
        for a, b in zip(direct, old):
            formula_err = max(formula_err, abs(a - b) / max(1.0, abs(b)))
        for which in (0, 3):
            for axis in range(4):
                h = 2e-6
                lo, hi = list(point), list(point)
                lo[axis] -= h
                hi[axis] += h
                fd = (selectors(tuple(hi))[0][0 if which == 0 else 1]
                      - selectors(tuple(lo))[0][0 if which == 0 else 1]) / (2*h)
                got = midpoint(ad[which].g[axis])
                deriv_err = max(deriv_err,
                                abs(fd - got) / max(1.0, abs(fd), abs(got)))
        # Same physical angles, now written as u=s-beta.
        relative_point = (point[0], point[1], point[2] - point[0], point[3])
        relative = selectors_mode(relative_point, True)[0]
        for a, b in zip(direct, relative):
            chart_err = max(chart_err, abs(a - b) / max(1.0, abs(a), abs(b)))
        for mode, q in ((False, point), (True, relative_point)):
            row, adrow = rows_mode(q, mode)
            for which in (0, 1):
                for axis in range(4):
                    h = 2e-6
                    lo, hi = list(q), list(q)
                    lo[axis] -= h
                    hi[axis] += h
                    fd = (rows_mode(tuple(hi), mode)[0][which]
                          - rows_mode(tuple(lo), mode)[0][which]) / (2 * h)
                    got = midpoint(adrow[which].g[axis])
                    row_deriv_err = max(
                        row_deriv_err,
                        abs(fd - got) / max(1.0, abs(fd), abs(got)))

    for point in samples[:50]:
        radii = (2e-5, 2e-5, 3e-5, 3e-5)
        box = tuple(v for pair in zip(
            (point[i] - radii[i] for i in range(4)),
            (point[i] + radii[i] for i in range(4))) for v in pair)
        mv = D.selector_mean_value({"coord": "absolute"}, *box)
        if mv is None:
            continue
        for k in range(9):
            q = tuple(point[i] + radii[i] * (2 * ((k * (i + 3)) % 11) / 10 - 1)
                      for i in range(4))
            val = selectors(q)[0]
            for got, enclosure in zip(val, mv):
                if not (J.lo(enclosure) <= got <= J.hi(enclosure)):
                    contain_fail += 1

    # Coordinate-row mean-value and Taylor-collar containment.
    for point in samples[:40]:
        for mode in (False, True):
            q0 = point if not mode else (
                point[0], point[1], point[2] - point[0], point[3])
            radii = (1e-5, 1e-5, 2e-5, 2e-5)
            box = tuple(v for pair in zip(
                (q0[i] - radii[i] for i in range(4)),
                (q0[i] + radii[i] for i in range(4))) for v in pair)
            mv = D.rowform_mean_value({"coord": "relative" if mode
                                       else "absolute"}, *box)
            for k in range(5):
                q = tuple(q0[i] + radii[i] *
                          (2 * ((k * (i + 2) + i) % 9) / 8 - 1)
                          for i in range(4))
                vals = rows_mode(q, mode)[0]
                if mv is not None:
                    for got, enclosure in zip(vals, mv):
                        if not contains(enclosure, got):
                            row_contain_fail += 1

            image = D.coordinate_krawczyk_image(
                {"coord": "relative" if mode else "absolute"}, *box)
            # The center need not be a zero, so only falsify a claimed empty
            # image when a numerical root is not available.  A contraction
            # must at least remain a subbox structurally.
            if image is not None and image[0] == "contract":
                z = image[1]
                if not (box[4] <= z[4] <= z[5] <= box[5]
                        and box[6] <= z[6] <= z[7] <= box[7]):
                    kraw_contraction_fail += 1

    for k in range(80):
        b = rng.uniform(0.15, 3.0)
        y = rng.uniform(-0.95, 0.9)
        s = rng.uniform(0.2, 6.0)
        gap = rng.uniform(1e-4, 0.7)
        for mode in (False, True):
            u = s - b if mode else s
            radii = (2e-5, 2e-5, 3e-5, min(2e-5, gap / 4))
            center = (b, y, u, gap)
            box = tuple(v for pair in zip(
                (center[i] - radii[i] for i in range(4)),
                (center[i] + radii[i] for i in range(4))) for v in pair)
            spec = {"coord": "relative" if mode else "absolute"}
            collar = D.collar_rowform(
                spec, D.mk(*box[0:2]), D.mk(*box[2:4]),
                D.mk(*box[4:6]), D.mk(*box[6:8]))
            exact = rows_mode(center, mode)[0]
            for got, enclosure in zip(exact, collar):
                if not contains(enclosure, got / gap ** 3):
                    collar_contain_fail += 1
        G = X.thin(gap)
        a, bq, q = D.small_gap_abq(G)
        with X.mpmath.workdps(100):
            # X.thin(float) encloses that exact binary float, not its shorter
            # decimal rendering; evaluate the reference at the same number.
            g = X.mpmath.mpf(gap)
            phi = (X.mpmath.pi - g) * X.mpmath.cos(g) + X.mpmath.sin(g)
            h = (X.mpmath.pi - g) * X.mpmath.sin(g)
            exact_a = (X.mpmath.pi - phi) / (g * g)
            exact_b = h / g
            exact_q = (2 * exact_a - exact_b) / g
        for enclosure, exact in ((a, exact_a), (bq, exact_b), (q, exact_q)):
            if not contains(enclosure, exact):
                quotient_contain_fail += 1

    print("samples=%d formula_relerr=%.3e derivative_relerr=%.3e "
          "chart_relerr=%.3e row_derivative_relerr=%.3e" %
          (len(samples), formula_err, deriv_err, chart_err, row_deriv_err))
    print("literal_relerr=%.3e selector_mv=%d row_mv=%d kraw=%d collar=%d "
          "quotients=%d failures" %
          (literal_err, contain_fail, row_contain_fail,
           kraw_contraction_fail, collar_contain_fail,
           quotient_contain_fail))
    if (formula_err > 2e-11 or deriv_err > 2e-5 or chart_err > 2e-11
            or row_deriv_err > 2e-5 or literal_err > 2e-11 or contain_fail
            or row_contain_fail or kraw_contraction_fail or collar_contain_fail
            or quotient_contain_fail):
        print("GATE FAILED")
        return 1
    print("GATE PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
