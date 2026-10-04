"""
Float self-check of the PROOF CHAIN of the separatedNumJ strict sign law
(Lean: Planar/SignedN2/FreeMassJ/SeparatedCount/SignLaw*.lean, 2026-08-08).

This is the numerical validation that preceded the transcription (validity
record, formalization-conventions section 1v): every identity and every
inequality the Lean proof uses, checked on dense random grids, PLUS the
final sign law re-derived on an independent grid.  Validity domain: float64,
all sample points kept >= 1e-3 away from the degenerate boundary lattice
(at the boundary itself the quantities vanish to order up to 4 and float64
cancellation manufactures spurious signs; the small-t behaviour of S, N, R2
was separately confirmed at 40 digits with mpmath during the campaign).

Checks (matching the Lean decomposition):
  1  termwise 1-D weights R, R2 >= 0; S > 0; norm defect N > 0 on (0, pi)
  2  Wronskian branch identities (x >= D and x <= D) vs direct W'h - Wh'
  3  Wronskian positivity on (0, pi)^2
  4  R2/sin^2 cross monotonicity (the g-mono chain input)
  5  fold split num = hJ(D) T+ + hJ(y) N  (y = 2pi - x)
  6  T+ segment-core identities (below/above the fold)
  7  T+ >= endpoint bounds (0 below the fold, pi*S(s-pi) above)
  8  odd reflection num(2pi-D, 2pi-x) = -num(D, x)
  9  the sign law itself on 5*10^5 random points
"""
import numpy as np

pi = np.pi
rng = np.random.default_rng(0)


def wrap(t):
    return np.arccos(np.cos(t))


def K(t):
    a = wrap(t)
    return (pi - a) * np.cos(t) + np.sin(a)


def h(t):
    a = wrap(t)
    return np.sin(t) * (pi - a)


def num(D, x):
    return (pi * K(x - D) - K(D) * K(x)) * h(D) - (pi ** 2 - K(D) ** 2) * h(x)


S = lambda t: np.sin(t) - t * np.cos(t)
St = lambda D: np.sin(D) + (2 * pi - D) * np.cos(D)
R = lambda x: (pi - x) ** 2 - np.sin(x) ** 2 + pi * ((pi - x) - np.sin(x) * np.cos(x))
R2 = lambda x: pi * (x - np.sin(x) * np.cos(x)) - (x ** 2 - np.sin(x) ** 2)
N = lambda D: pi ** 2 - K(D) ** 2 - h(D) ** 2
W = lambda x, D: pi * K(x - D) - K(D) * K(x)


def main():
    x = np.linspace(1e-3, pi - 1e-3, 300000)
    assert (R(x) >= 0).all() and (R2(x) >= 0).all()
    assert (S(x) > 0).all() and (N(x) > 0).all()

    D = rng.uniform(1e-3, pi - 1e-3, 200000)
    X = rng.uniform(1e-3, pi - 1e-3, 200000)
    Wr_formula = np.where(
        X >= D,
        S(D) * R(X) + pi * D * np.sin(D) * np.sin(X) ** 2,
        -St(D) * R2(X) + pi * (2 * pi - D) * np.sin(D) * np.sin(X) ** 2)
    Wx = -pi * h(X - D) + K(D) * h(X)
    hp = K(X) - 2 * np.sin(X)
    Wr_direct = Wx * h(X) - W(X, D) * hp
    assert np.allclose(Wr_formula, Wr_direct, atol=1e-9)
    assert (Wr_formula > 0).all()

    mask = X < D
    gx = R2(X) / np.sin(X) ** 2
    gD = R2(D) / np.sin(D) ** 2
    assert ((gx < gD) | ~mask).all()

    D2 = rng.uniform(1e-3, pi - 1e-3, 200000)
    X2 = rng.uniform(pi + 1e-3, 2 * pi - 1e-3, 200000)
    Y = 2 * pi - X2
    Tplus = pi * K(Y + D2) - K(D2) * K(Y) + h(D2) * h(Y)
    assert np.allclose(num(D2, X2), h(D2) * Tplus + h(Y) * N(D2), atol=1e-8)
    sE = D2 + Y
    E = lambda t, s: (t * (s - t) * np.sin(t) * np.sin(s - t)
                      - S(t) * S(s - t))
    Tp_id = np.where(sE <= pi, E(D2, sE), E(D2, sE) + 2 * pi * S(sE - pi))
    assert np.allclose(Tplus, Tp_id, atol=1e-9)
    m1 = sE <= pi
    assert (Tp_id[m1] > 0).all()
    assert (Tp_id[~m1] >= pi * S(sE[~m1] - pi) - 1e-9).all()

    D3 = rng.uniform(1e-3, 2 * pi - 1e-3, 200000)
    X3 = rng.uniform(1e-3, 2 * pi - 1e-3, 200000)
    assert np.allclose(num(2 * pi - D3, 2 * pi - X3), -num(D3, X3), atol=1e-8)

    Dg = rng.uniform(1e-4, 2 * pi - 1e-4, 500000)
    Xg = rng.uniform(1e-4, 2 * pi - 1e-4, 500000)
    v = num(Dg, Xg)
    pred = np.where(Xg > Dg, 1, -1)
    bad = (np.sign(v) != pred) & (np.abs(Xg - Dg) > 1e-3)
    assert bad.sum() == 0
    print("ALL SIGN-LAW CHAIN CHECKS PASSED")


if __name__ == "__main__":
    main()
