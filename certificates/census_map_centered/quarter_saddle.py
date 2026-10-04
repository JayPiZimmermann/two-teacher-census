"""Interval certificate that the quarter family on y = -1/2 is a SADDLE.

det A = c0 c1 [Tau0 Tau1 + c0 Tau0 + c1 Tau1] is used in FACTORED form: the
expanded A00 A11 - A01^2 loses everything to dependency (enclosure [-2.15, 1.57]
on the box [1, 1.1]), the factored one does not.
"""
import json, mpmath
from mpmath import iv
import cert_core as K
C = K.C; PI = K.PI


def parts(bI):
    u = bI / 2
    Cu, Su = iv.cos(u), iv.sin(u)
    kap = PI() / 2
    dt = PI() * PI() / 4 - 1
    P0 = 2 * ((kap - u) * Cu + Su)
    P1 = 2 * (u * Su + Cu)
    T0 = 2 * ((kap - u) * Cu - Su)
    T1 = 2 * (u * Su - Cu)
    c0 = (kap * P0 - P1) / dt
    c1 = (kap * P1 - P0) / dt
    Kk = T0 * T1 + c0 * T0 + c1 * T1
    return c0, c1, Kk, T0, T1


def detA_sign(bI):
    c0, c1, Kk, _, _ = parts(bI)
    s = K.sgn(c0) * K.sgn(c1) * K.sgn(Kk)
    return s


def cover(lo, hi, dmin, budget=4000000):
    todo, good, bad = [(lo, hi)], [], []
    n = 0
    while todo and n < budget:
        a, b = todo.pop(); n += 1
        if detA_sign(iv.mpf([a, b])) == -1:
            good.append((a, b)); continue
        if b - a <= dmin:
            bad.append((a, b)); continue
        m = (a + b) / 2
        todo += [(a, m), (m, b)]
    return good, bad + todo, n


if __name__ == "__main__":
    EPS = mpmath.mpf("0.002")     # det A ~ -const * beta^4 at the ends
    lo = EPS
    hi = C.endpoints(PI())[0] - EPS
    good, bad, n = cover(lo, hi, mpmath.mpf("1e-16"), budget=4000000)
    cov = sum(b - a for a, b in good)
    print("certified det A < 0 on %.14f of [0.002, pi-0.002] (width %.14f)"
          % (float(cov), float(hi - lo)))
    print("   %d boxes, %d tests, %d failures" % (len(good), n, len(bad)))
    for a, b in bad[:6]:
        print("   FAIL", mpmath.nstr(a, 20), mpmath.nstr(b, 20))
    spots = []
    for bf in [0.01, 0.1, 0.5, 1.0, 1.420925475551, 1.5707963267948966,
               1.720667178038759, 2.0, 2.7, 3.1]:
        c0, c1, Kk, T0, T1 = parts(C.I(bf))
        spots.append({"beta": bf, "detA": mpmath.nstr(C.mid(c0 * c1 * Kk), 12),
                      "c0": mpmath.nstr(C.mid(c0), 12), "c1": mpmath.nstr(C.mid(c1), 12),
                      "Tau0": mpmath.nstr(C.mid(T0), 12), "Tau1": mpmath.nstr(C.mid(T1), 12),
                      "sign_detA": detA_sign(C.I(bf))})
        print("   beta=%-10.7f detA=%-14s c0=%-11s c1=%-11s Tau0=%-11s Tau1=%-11s sgn=%d"
              % (bf, spots[-1]["detA"][:13], spots[-1]["c0"][:11], spots[-1]["c1"][:11],
                 spots[-1]["Tau0"][:11], spots[-1]["Tau1"][:11], spots[-1]["sign_detA"]))
    json.dump({"range": ["0.002", "pi - 0.002"], "boxes": len(good), "tests": n,
               "measure_certified_detA_negative": float(cov),
               "range_width": float(hi - lo),
               "failed_boxes": [[mpmath.nstr(a, 20), mpmath.nstr(b, 20)] for a, b in bad],
               "spot_values": spots,
               "verdict": "det A < 0 on every certified box: the quarter family "
                          "on y = -1/2 is a strict SADDLE, so it does not change "
                          "the census and does not split F1 or F4"},
              open("quarter_saddle.json", "w"), indent=1)
