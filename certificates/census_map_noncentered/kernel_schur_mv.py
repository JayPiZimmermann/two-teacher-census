"""Mean-value (first-order) enclosure of the mass-free Schur determinant.

WHY.  Measured across the 45 non-exact-fit certified families, the LITERAL
enclosure of `generalJKernelTeacherSchurDetT` has a width growing linearly in
the box side (9.5e7 * side at the F4 family, against a value of -1541), so the
sign test fires only below a family-dependent side: median 2e-4, with a tail
of five families unresolved down to 1e-7.  That tail is the whole cost of the
row-7a sweep.  The cost driver is the enclosure's DEPENDENCY SLOPE, not the
margin -- every one of those families has a relative margin of at least 0.205
-- and a first-order form is exactly what attacks a slope.

HOW.  The determinant is a long composition (num -> kernel vector -> V, L, N0,
N1 -> T00, T11, T01 -> det), so its derivative is taken by FORWARD-MODE
AUTOMATIC DIFFERENTIATION over intervals rather than by hand: each quantity is
carried as a triple `(value, d/ds, d/dD)` in the census coordinates
`th0 = s + D`, `th1 = s`, and the kernel atoms seed it with the tree's own
derivative identities `phi' = -h`, `h' = sA`, and `d|sin t|/dt = cos p` for
the branch coordinate `p` in `[0, pi]`.  Hand-deriving an eight-deep chain is
where a sign error hides, and this certificate has already paid once for an
invalid derivative enclosure.

GATE: `kernel_schur_mv_gate.py` checks the AD derivative against central
differences and checks that the mean-value form encloses the literal one at
sample points.  Nothing here is used until that passes.
"""
import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS

PIv = KS.PIv


class D2(object):
    """(value, d/ds, d/dD) over intervals; forward-mode AD."""
    __slots__ = ("v", "s", "d")

    def __init__(self, v, s=None, d=None):
        z = v * 0
        self.v = v
        self.s = z if s is None else s
        self.d = z if d is None else d

    def __add__(self, o):
        o = _lift(o, self.v)
        return D2(self.v + o.v, self.s + o.s, self.d + o.d)

    __radd__ = __add__

    def __neg__(self):
        return D2(-self.v, -self.s, -self.d)

    def __sub__(self, o):
        return self + (-_lift(o, self.v))

    def __rsub__(self, o):
        return _lift(o, self.v) + (-self)

    def __mul__(self, o):
        o = _lift(o, self.v)
        return D2(self.v * o.v,
                  self.s * o.v + self.v * o.s,
                  self.d * o.v + self.v * o.d)

    __rmul__ = __mul__


def _lift(o, proto):
    return o if isinstance(o, D2) else D2(proto * 0 + o)


def _atom_D2(T, ds, dd):
    """phi, h, |sin| at an angle whose (d/ds, d/dD) are (ds, dd)."""
    A = KS._atoms(T)
    absin = KS._abs_sin(T)
    dabs = _dabs_sin(T)
    phi = D2(A.phi, (-A.h) * ds, (-A.h) * dd)
    h = D2(A.h, A.sA * ds, A.sA * dd)
    a = D2(absin, dabs * ds, dabs * dd)
    return phi, h, a


def _dabs_sin(T):
    """`d|sin t|/dt`, with the tree's branch convention.

    On branch `n`, `|sin t| = e*sin(p)` with `e=+1` on even branches and
    `e=-1` on odd branches.  Since `dp/dt=-1`, its derivative is
    `-e*cos(p)`.  Omitting `e` makes both the literal selector and this AD form
    wrong on every odd branch while leaving a self-comparison gate falsely
    green."""
    out = None
    for piece, n in J.branch_pieces(T):
        e = 1 if n % 2 == 0 else -1
        v = -(e * X.iv.cos(J.p_of(piece, n)))
        out = v if out is None else J.hull(out, v)
    return out


def detT_D2(beta, T0, T1, Dm):
    """(det, d det/ds, d det/dD) with `beta` thin and `th0 = s + D`."""
    one = Dm * 0 + 1
    zero = Dm * 0
    phiD, hD, aD = _atom_D2(Dm, zero, one)
    phi0, h0, a0 = _atom_D2(T0, one, one)
    phi1, h1, a1 = _atom_D2(T1, one, zero)
    phi0b, h0b, a0b = _atom_D2(T0 - beta, one, one)
    phi1b, h1b, a1b = _atom_D2(T1 - beta, one, zero)
    # NOTE: mpmath's interval type RAISES on an unknown right operand instead
    # of returning NotImplemented, so `__rmul__` never fires and a bare
    # `PIv * <D2>` dies.  Every constant is lifted to a D2 (with zero
    # derivative) before it meets one.
    PI = D2(PIv)
    M = D2(PIv * PIv) - phiD * phiD
    # num(D, th0) uses phi(th0 - D) = phi(th1); num(D, th0-beta) uses phi(th1-beta)
    S1 = (PI * phi1 - phiD * phi0) * hD - M * h0
    S0 = -((PI * phi1b - phiD * phi0b) * hD - M * h0b)
    V0 = S0 * phi0 + S1 * phi0b
    V1 = S0 * phi1 + S1 * phi1b
    L0 = S0 * a0 + S1 * a0b
    L1 = S0 * a1 + S1 * a1b
    N0 = PI * V0 - phiD * V1
    N1 = PI * V1 - phiD * V0
    K = 2 * aD - phiD
    A = M * N0 * N1 * K
    P = M * M * N0 * (2 * L0 - V0) + PI * hD * hD * N0 * N0
    Q = M * M * N1 * (2 * L1 - V1) + PI * hD * hD * N1 * N1
    R = phiD * hD * hD * N0 * N1
    det = -(A * (P + Q + 2 * R)) + P * Q - R * R
    return det.v, det.s, det.d


def detT_meanvalue(beta, s_lo, s_hi, d_lo, d_hi):
    """Mean-value enclosure of the determinant over the (s, D) box."""
    T1 = X.iv.mpf([s_lo, s_hi])
    Dm = X.iv.mpf([d_lo, d_hi])
    _, gs, gd = detT_D2(beta, T1 + Dm, T1, Dm)
    sm = X.thin(J.mid(T1))
    dm = X.thin(J.mid(Dm))
    c, _, _ = detT_D2(beta, sm + dm, sm, dm)
    return c + gs * (T1 - sm) + gd * (Dm - dm)
