"""
The QUARTER locus: a fourth mass-free determinant the certificate does not cover,
and the proof that it does not split a face.

A separated critical point with D = th0 - th1 = pi/2 (mod pi) has H(D) = 0, so
the torque rows cannot be solved for the masses; instead BOTH students must sit
at zeros of the TEACHER torque, and the masses come from the radial rows.  Such a
family exists exactly where two teacher torque roots are a quarter turn apart,
i.e. on the zero set of

    G(beta, t) = H(t-beta) H(t+pi/2) - H(t) H(t+pi/2-beta).

G is a mass-free determinant like Wtau, Wpot, Wwgt, but it is NOT among the three
of the certificate, and the shipped classification lists the family it carries
(`centeredQuarterRows`) only when T.s0 === T.s1 EXACTLY -- a float64 equality
that `massesAt(-0.5)` never satisfies, so the shipped map never draws it.

Part 1 enumerates the zeros of G over beta and shows the locus is exactly
    {y = 0}  u  {y = +-1}  u  {y = -1/2},
the first two already being strata of the arrangement and the third being new.

Part 2 certifies, in interval arithmetic, that the family on y = -1/2 is a
SADDLE for every beta in the covered range, so it changes the family inventory
but not the CENSUS, and therefore does not split the faces F1 / F4 that the line
y = -1/2 runs through.

Closed forms on y = -1/2 (s0 = s1; the type is invariant under positive scaling
of the teacher, so take s0 = s1 = 1).  With u = beta/2, C = cos u, S = sin u:
    th0 = u,  th1 = u + pi/2,  D = -pi/2,  phi(D) = 1,  H(D) = 0,  H'(D) = -1
    P0 = 2[(pi/2-u)C + S]      P1 = 2[uS + C]
    Tau0 = 2[(pi/2-u)C - S]    Tau1 = 2[uS - C]
    c0 = (kap P0 - P1)/dt,  c1 = (kap P1 - P0)/dt,  dt = pi^2/4 - 1,  kap = pi/2
Because H(D) = 0 the mass/angle coupling block B vanishes and the Hessian is
block diagonal: Gram (PD) plus
    A = [[c0 Tau0 + c0 c1, -c0 c1], [-c0 c1, c1 Tau1 + c0 c1]],
so the point is a local minimum iff A is positive semidefinite, and
    det A = c0 c1 [ Tau0 Tau1 + c0 Tau0 + c1 Tau1 ].
Note Tau0 = 0 exactly at beta = beta1* and Tau1 = 0 exactly at beta = beta2*:
the two lens endpoints are where this family's angular Hessian loses a diagonal
entry.
"""
import json
import mpmath
from mpmath import iv
import cert_core as K

C = K.C
PI = K.PI


def quarter_A(betaI):
    """(detA, A00, c0, c1) enclosures on y = -1/2 for beta in betaI."""
    u = betaI / 2
    Cu, Su = iv.cos(u), iv.sin(u)
    kap = PI() / 2
    dt = PI() * PI() / 4 - 1
    P0 = 2 * ((kap - u) * Cu + Su)
    P1 = 2 * (u * Su + Cu)
    T0 = 2 * ((kap - u) * Cu - Su)
    T1 = 2 * (u * Su - Cu)
    c0 = (kap * P0 - P1) / dt
    c1 = (kap * P1 - P0) / dt
    A00 = c0 * T0 + c0 * c1
    A11 = c1 * T1 + c0 * c1
    A01 = -c0 * c1
    detA = A00 * A11 - A01 * A01
    return detA, A00, c0, c1, T0, T1


def cover_saddle(lo, hi, dmin):
    todo = [(lo, hi)]
    good, bad = [], []
    while todo:
        a, b = todo.pop()
        BI = iv.mpf([a, b])
        detA = quarter_A(BI)[0]
        if K.sgn(detA) == -1:
            good.append((a, b))
            continue
        if b - a <= dmin:
            bad.append((a, b))
            continue
        m = (a + b) / 2
        todo += [(a, m), (m, b)]
    return good, bad


# ---------------------------------------------------------------- part 1

def H(t):
    x = mpmath.fmod(t, mpmath.pi)
    if x < 0:
        x += mpmath.pi
    return (mpmath.pi / 2 - x) * mpmath.sin(x)


def G(beta, t):
    return H(t - beta) * H(t + mpmath.pi / 2) - H(t) * H(t + mpmath.pi / 2 - beta)


def ratio_y(beta, t):
    s0, s1 = -H(t - beta), H(t)
    a = mpmath.atan2(s0, s1)
    a = mpmath.fmod(a, mpmath.pi)
    if a <= 0:
        a += mpmath.pi
    return 2 * a / mpmath.pi - 1


def enumerate_locus(nbeta=120, N=6000):
    mpmath.mp.dps = 40
    ys = {}
    for i in range(1, nbeta):
        beta = mpmath.pi * mpmath.mpf(i) / nbeta
        prev_t, prev = mpmath.mpf(0), G(beta, mpmath.mpf(0))
        found = []
        for j in range(1, N + 1):
            t = mpmath.pi * mpmath.mpf(j) / N
            v = G(beta, t)
            if prev * v < 0:
                a, b = prev_t, t
                for _ in range(120):
                    m = (a + b) / 2
                    if G(beta, a) * G(beta, m) <= 0:
                        b = m
                    else:
                        a = m
                found.append((a + b) / 2)
            prev_t, prev = t, v
        found += [mpmath.mpf(0), mpmath.pi / 2]      # the two endpoint zeros
        for t in found:
            yv = float(ratio_y(beta, t))
            key = round(yv, 9)
            ys.setdefault(key, 0)
            ys[key] += 1
    return ys


def main():
    rep = {}
    print("part 1: the y-values of the quarter locus over 119 betas")
    ys = enumerate_locus()
    rep["locus_y_values"] = {str(k): v for k, v in sorted(ys.items())}
    for k, v in sorted(ys.items()):
        print("   y = %+.9f   hit %d times" % (k, v))

    print("\npart 2: certified sign of det A on y = -1/2")
    dmin = mpmath.mpf("1e-12")
    good, bad = cover_saddle(mpmath.mpf("1e-6"),
                             C.endpoints(PI())[0] - mpmath.mpf("1e-6"), dmin)
    covered = sum(b - a for a, b in good)
    rep["saddle_certificate"] = {
        "range": ["1e-6", "pi - 1e-6"],
        "boxes": len(good),
        "measure_certified_detA_negative": float(covered),
        "failed_boxes": [[mpmath.nstr(a, 20), mpmath.nstr(b, 20)] for a, b in bad],
        "verdict": "det A < 0 (indefinite angular block) on every certified box, "
                   "so the quarter family is a strict SADDLE there and adds "
                   "nothing to the census"}
    print("   certified det A < 0 on %.9f of (1e-6, pi-1e-6), %d boxes, %d failures"
          % (float(covered), len(good), len(bad)))
    for a, b in bad[:10]:
        print("      FAILED", mpmath.nstr(a, 20), mpmath.nstr(b, 20))

    # spot values, for the record
    spots = []
    for bf in [0.01, 0.1, 0.5, 1.0, 1.420925475551, 1.5707963267948966,
               1.720667178038759, 2.0, 2.7, 3.1]:
        detA, A00, c0, c1, T0, T1 = quarter_A(C.I(bf))
        spots.append({"beta": bf,
                      "detA": mpmath.nstr(C.mid(detA), 12),
                      "A00": mpmath.nstr(C.mid(A00), 12),
                      "c0": mpmath.nstr(C.mid(c0), 12),
                      "c1": mpmath.nstr(C.mid(c1), 12),
                      "Tau0": mpmath.nstr(C.mid(T0), 12),
                      "Tau1": mpmath.nstr(C.mid(T1), 12),
                      "sign_detA": K.sgn(detA)})
    rep["spot_values"] = spots
    print("\n   beta        detA          c0        c1        Tau0      Tau1")
    for s in spots:
        print("   %-11.7f %-13s %-9s %-9s %-9s %-9s" %
              (s["beta"], s["detA"][:12], s["c0"][:9], s["c1"][:9],
               s["Tau0"][:9], s["Tau1"][:9]))

    with open("quarter.json", "w") as fh:
        json.dump(rep, fh, indent=1)


if __name__ == "__main__":
    main()
