"""Adaptive interval DFS for the NONCENTERED SELECTED minimum map.

This is the numerical consumer of Lean's
``CountConstancyJ/MinimumSelectorJ`` theorem.  It sweeps a named box in

    (beta, y, s, D),  theta1 = s, theta0 = s + D,

or, with ``coord=relative``, in the teacher-1-aligned chart

    (beta, y, u, D),  theta1 = beta + u, theta0 = beta + u + D.

and discharges every leaf by one of the following one-sided facts:

``MAP_PLAIN`` / ``MAP_CENTERED``
    A component of the eliminated angle map misses zero, so there is no root.
``MAP_AUX``
    Either a parameter-uniform Krawczyk image misses the angle box, or the
    separated sign law excludes it after both Cramer numerators are certified
    nonzero.
``MAP_COLLAR``
    A Taylor-normalized angle equation misses zero near student coincidence.
``NEG_T00``
    The Cramer-cleared first Schur pivot is strictly negative, so no root in
    the box belongs to the selected set.
``NEG_DET``
    The Cramer-cleared Schur determinant is strictly negative, so no root in
    the box belongs to the selected set.
``POS_DET``
    The Cramer-cleared determinant is strictly positive.  Any selected root
    in the box is therefore regular, exactly as required by
    ``CensusMinimumSelectedRegularOnJ``.

Nothing is asked of unselected roots.  In particular a saddle-saddle fold is
closed by ``NEG_T00`` or ``NEG_DET`` and is not a wall of this certificate.
Only a zero with nonnegative first pivot and determinant zero can survive all
four alternatives.

Certificate format: a preorder bisection forest.  Root cells are a specified
equal partition of the full-period third coordinate.  At a node, ``1`` means
internal and its two children follow; ``0`` means leaf and is followed by a
three-bit verdict.  The split is deterministic: isolate an interval endpoint
of ``u=n*pi`` in the relative chart when one is present, otherwise bisect the
widest of the four coordinates; child zero is the lower piece.  Bits are
packed MSB first.

The search is not the claim.  ``minimum_dfs_check.py`` regenerates every root
and child and recomputes every verdict.

Usage:
  python3 minimum_dfs.py b0 b1 y0 y1 delta dmax tag \
      n_parts lo_idx hi_idx [budget] [minw] [absolute|relative]
"""
import base64
import json
import os
import sys
import time

import mpmath

import census_cert as X
import face_bb_signlaw as F
import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))

VERDICTS = (
    "MAP_PLAIN",
    "MAP_CENTERED",
    "MAP_AUX",
    "MAP_COLLAR",
    "NEG_T00",
    "NEG_DET",
    "POS_DET",
    "UNDECIDED",
)
CODE = {name: i for i, name in enumerate(VERDICTS)}
NAME = {i: name for name, i in CODE.items()}


def mk(lo, hi):
    def scalar(z):
        return z if isinstance(z, mpmath.mpf) else mpmath.mpf(repr(z))
    return X.iv.mpf([scalar(lo), scalar(hi)])


def root_box(spec):
    """The full region, regenerated from the spec alone."""
    two_pi = float(J.hi(2 * X.PIv))
    return (
        spec["b0"], spec["b1"], spec["y0"], spec["y1"],
        spec["seam"], spec["seam"] + two_pi,
        spec["delta"], spec["dmax"],
    )


def root_cells(spec):
    """Equal ``s`` partition of the root box."""
    bl, bh, yl, yh, sl, sh, dl, dh = root_box(spec)
    n = spec["n_parts"]
    return [
        (bl, bh, yl, yh,
         sl + (sh - sl) * i / n,
         sl + (sh - sl) * (i + 1) / n,
         dl, dh)
        for i in range(n)
    ]


def partition_covers(spec):
    root = root_box(spec)
    cells = root_cells(spec)
    if not cells or cells[0][4] != root[4] or cells[-1][5] != root[5]:
        return False
    for i, cell in enumerate(cells):
        if cell[:4] != root[:4] or cell[6:] != root[6:]:
            return False
        if i and cells[i - 1][5] != cell[4]:
            return False
    return True


def split_location(spec, bx):
    """Deterministic branch-aligned cut, otherwise widest-axis midpoint."""
    # In the relative chart the teacher-1 atom is evaluated exactly at u.
    # Isolate every u=n*pi branch boundary using the outward endpoints of
    # mpmath's interval pi.  This prevents a single irrational kink from
    # forcing irrelevant beta/y subdivisions forever; the tiny [lo(pi),
    # hi(pi)] sliver remains explicitly covered.
    if spec.get("coord") == "relative":
        sl, sh = bx[4], bx[5]
        n0 = int(mpmath.floor(mpmath.mpf(sl) / mpmath.pi)) - 1
        n1 = int(mpmath.ceil(mpmath.mpf(sh) / mpmath.pi)) + 1
        for n in range(n0, n1 + 1):
            z = X.iv.mpf(n) * X.PIv
            for cut in J.endpoints(z):
                if sl < cut < sh:
                    return 2, cut
    widths = tuple(bx[2 * i + 1] - bx[2 * i] for i in range(4))
    axis = widths.index(max(widths))
    return axis, (bx[2 * axis] + bx[2 * axis + 1]) / 2


def child(spec, bx, bit):
    """Apply the prescribed split; zero is the lower child."""
    axis, mid = split_location(spec, bx)
    out = list(bx)
    out[2 * axis + (1 if bit == 0 else 0)] = mid
    return tuple(out)


def volume(bx):
    out = 1.0
    for i in range(4):
        out *= bx[2 * i + 1] - bx[2 * i]
    return out


class D4:
    """Interval value and gradient in ``(beta, y, s, D)``."""
    __slots__ = ("v", "g")

    def __init__(self, value, grad=None):
        self.v = value
        zero = value * 0
        self.g = (zero, zero, zero, zero) if grad is None else tuple(grad)

    def __add__(self, other):
        other = lift(other, self.v)
        return D4(self.v + other.v,
                  tuple(self.g[i] + other.g[i] for i in range(4)))

    __radd__ = __add__

    def __neg__(self):
        return D4(-self.v, tuple(-z for z in self.g))

    def __sub__(self, other):
        return self + (-lift(other, self.v))

    def __rsub__(self, other):
        return lift(other, self.v) + (-self)

    def __mul__(self, other):
        other = lift(other, self.v)
        return D4(
            self.v * other.v,
            tuple(self.g[i] * other.v + self.v * other.g[i]
                  for i in range(4)),
        )

    __rmul__ = __mul__


def lift(value, proto):
    return value if isinstance(value, D4) else D4(proto * 0 + value)


def variable(value, axis):
    zero = value * 0
    grad = [zero, zero, zero, zero]
    grad[axis] = zero + 1
    return D4(value, grad)


def atom(a):
    """The kernel atoms composed with an interval AD angle."""
    q, smooth = J.atoms(a.v)
    return (
        D4(q.phi, tuple(q.dphi * z for z in a.g)),
        D4(q.h, tuple(q.dh * z for z in a.g)),
        D4(q.asin, tuple(q.dasin * z for z in a.g)),
        smooth,
    )


def census_masses(y):
    psi = (y + 1) * D4(X.PIv) * mpmath.mpf("0.5")
    sinv = X.iv.sin(psi.v)
    cosv = X.iv.cos(psi.v)
    return (
        D4(sinv, tuple(cosv * z for z in psi.g)),
        D4(cosv, tuple(-sinv * z for z in psi.g)),
    )


def cramer_selectors_d4(Bv, Yv, Sv, Dv, relative=False):
    """Lean's four Cramer-cleared Schur expressions, with interval AD.

    Returns ``(T00, T11, T01, det, smooth)``.  No student-weight or radial
    determinant division occurs.
    """
    B = variable(Bv, 0)
    Y = variable(Yv, 1)
    S = variable(Sv, 2)
    D = variable(Dv, 3)
    T1 = B + S if relative else S
    T0 = T1 + D
    U1 = S if relative else T1 - B
    U0 = S + D if relative else T0 - B
    s0, s1 = census_masses(Y)
    phiD, hD, absD, smoothD = atom(D)
    phi0, _, abs0, smooth0 = atom(T0)
    phi1, _, abs1, smooth1 = atom(T1)
    phi0b, _, abs0b, smooth0b = atom(U0)
    phi1b, _, abs1b, smooth1b = atom(U1)
    PI = D4(X.PIv)
    M = PI * PI - phiD * phiD
    P0 = s0 * phi0 + s1 * phi0b
    P1 = s0 * phi1 + s1 * phi1b
    L0 = s0 * abs0 + s1 * abs0b
    L1 = s0 * abs1 + s1 * abs1b
    N0 = PI * P0 - phiD * P1
    N1 = PI * P1 - phiD * P0
    band = 2 * absD - phiD
    # Exact ring factorization of Lean's expanded definitions.  Besides
    # being substantially sharper under interval evaluation, this exposes
    # the unavoidable N0*N1 factor of the determinant.
    A0 = M * N1 * band - M * M * (2 * L0 - P0) \
        - PI * hD * hD * N0
    A1 = M * N0 * band - M * M * (2 * L1 - P1) \
        - PI * hD * hD * N1
    C = -(M * band + phiD * hD * hD)
    T00 = N0 * A0
    T11 = N1 * A1
    T01 = N0 * N1 * C
    det = N0 * N1 * (A0 * A1 - N0 * N1 * C * C)
    return T00, T11, T01, det, all(
        (smoothD, smooth0, smooth1, smooth0b, smooth1b))


def coordinate_rowform_d4(Bv, Yv, Sv, Dv, relative=False):
    """The two division-free angle rows with interval AD in this chart."""
    B = variable(Bv, 0)
    Y = variable(Yv, 1)
    S = variable(Sv, 2)
    D = variable(Dv, 3)
    s0, s1 = census_masses(Y)
    T1 = B + S if relative else S
    T0 = T1 + D
    U1 = S if relative else T1 - B
    U0 = S + D if relative else T0 - B

    def load(Z0, Z1):
        p0, h0, _, z0 = atom(Z0)
        p1, h1, _, z1 = atom(Z1)
        return s0 * p0 + s1 * p1, -(s0 * h0 + s1 * h1), z0 and z1

    P0, A0, ok0 = load(T0, U0)
    P1, A1, ok1 = load(T1, U1)
    phi, h, _, okd = atom(D)
    PI = D4(X.PIv)
    M = PI * PI - phi * phi
    E0 = (PI * P1 - phi * P0) * h + M * A0
    E1 = (PI * P0 - phi * P1) * h - M * A1
    # The row map uses only phi and H.  Those atoms are globally C1 and the
    # branch derivatives agree at every k*pi endpoint; J.atoms already hulls
    # all one-sided values.  The smooth flags also track |sin|, whose
    # derivative does jump, but |sin| does not occur here.  Hence the row AD
    # remains a valid derivative enclosure across a branch endpoint.
    return E0, E1, True


def rowform_mean_value(spec, bl, bh, yl, yh, sl, sh, dl, dh):
    """First-order enclosure of both row equations in either chart."""
    boxes = (mk(bl, bh), mk(yl, yh), mk(sl, sh), mk(dl, dh))
    relative = spec.get("coord") == "relative"
    full = coordinate_rowform_d4(*boxes, relative)
    if not full[2]:
        return None
    mids = tuple(X.thin(J.mid(z)) for z in boxes)
    center = coordinate_rowform_d4(*mids, relative)

    def form(which):
        q = center[which].v
        for i in range(4):
            q += full[which].g[i] * (boxes[i] - mids[i])
        return q

    return form(0), form(1)


def coordinate_krawczyk_image(spec, bl, bh, yl, yh, sl, sh, dl, dh):
    """Parameter-uniform Krawczyk image in the two angle coordinates.

    ``(beta,y)`` stay interval parameters.  The inverse is formed at the
    all-thin midpoint, while the residual at the angle midpoint retains the
    whole parameter box.  Thus an empty intersection excludes a zero for
    every parameter in the box.  Otherwise every zero lies in the returned
    intersection with the angle box; this contraction is useful even without
    a uniqueness claim.
    """
    boxes = (mk(bl, bh), mk(yl, yh), mk(sl, sh), mk(dl, dh))
    relative = spec.get("coord") == "relative"
    full = coordinate_rowform_d4(*boxes, relative)
    if not full[2]:
        return None
    mids = tuple(X.thin(J.mid(z)) for z in boxes)
    point = coordinate_rowform_d4(*mids, relative)
    a00, a01 = J.mid(point[0].g[2]), J.mid(point[0].g[3])
    a10, a11 = J.mid(point[1].g[2]), J.mid(point[1].g[3])
    det = a00 * a11 - a01 * a10
    if det == 0:
        return None
    C = (X.thin(a11 / det), X.thin(-a01 / det),
         X.thin(-a10 / det), X.thin(a00 / det))
    # Residual at the angle midpoint, uniformly over the parameter box.
    F0, F1, ok = coordinate_rowform_d4(
        boxes[0], boxes[1], mids[2], mids[3], relative)
    if not ok:
        return None
    j00, j01 = full[0].g[2], full[0].g[3]
    j10, j11 = full[1].g[2], full[1].g[3]
    q00 = 1 - (C[0] * j00 + C[1] * j10)
    q01 = -(C[0] * j01 + C[1] * j11)
    q10 = -(C[2] * j00 + C[3] * j10)
    q11 = 1 - (C[2] * j01 + C[3] * j11)
    rs, rd = boxes[2] - mids[2], boxes[3] - mids[3]
    K0 = mids[2] - (C[0] * F0.v + C[1] * F1.v) + q00 * rs + q01 * rd
    K1 = mids[3] - (C[2] * F0.v + C[3] * F1.v) + q10 * rs + q11 * rd
    ns = (max(J.lo(K0), J.lo(boxes[2])),
          min(J.hi(K0), J.hi(boxes[2])))
    nd = (max(J.lo(K1), J.lo(boxes[3])),
          min(J.hi(K1), J.hi(boxes[3])))
    if ns[0] > ns[1] or nd[0] > nd[1]:
        return "empty", None
    return "contract", (bl, bh, yl, yh, ns[0], ns[1], nd[0], nd[1])


def selector_mean_value(spec, bl, bh, yl, yh, sl, sh, dl, dh, full=None):
    """First-order enclosures of ``T00`` and ``det`` over one box.

    The form is used only when every atom interval lies in one smooth
    half-branch.  Boxes crossing a kink fall back to literal evaluation.
    """
    boxes = (mk(bl, bh), mk(yl, yh), mk(sl, sh), mk(dl, dh))
    if full is None:
        full = cramer_selectors_d4(*boxes, spec.get("coord") == "relative")
    if not full[4]:
        return None
    mids = tuple(X.thin(J.mid(z)) for z in boxes)
    center = cramer_selectors_d4(*mids, spec.get("coord") == "relative")

    def form(which):
        q = center[which].v
        for i in range(4):
            q += full[which].g[i] * (boxes[i] - mids[i])
        return q

    return form(0), form(3)


def cramer_selectors(B, Y, S, D, relative=False):
    """Literal interval form of the four selectors, without AD overhead."""
    T1 = B + S if relative else S
    T0 = T1 + D
    U1 = S if relative else T1 - B
    U0 = S + D if relative else T0 - B
    s0, s1 = X.masses_at(Y)
    atom_d, _ = J.atoms(D)
    atom0, _ = J.atoms(T0)
    atom1, _ = J.atoms(T1)
    atom0b, _ = J.atoms(U0)
    atom1b, _ = J.atoms(U1)
    phi, h, abs_d = atom_d.phi, atom_d.h, atom_d.asin
    mass_det = X.PIv * X.PIv - phi * phi
    p0 = s0 * atom0.phi + s1 * atom0b.phi
    p1 = s0 * atom1.phi + s1 * atom1b.phi
    load0 = s0 * atom0.asin + s1 * atom0b.asin
    load1 = s0 * atom1.asin + s1 * atom1b.asin
    num0 = X.PIv * p0 - phi * p1
    num1 = X.PIv * p1 - phi * p0
    band = 2 * abs_d - phi
    a0 = mass_det * num1 * band \
        - mass_det * mass_det * (2 * load0 - p0) \
        - X.PIv * h * h * num0
    a1 = mass_det * num0 * band \
        - mass_det * mass_det * (2 * load1 - p1) \
        - X.PIv * h * h * num1
    c = -(mass_det * band + phi * h * h)
    t00 = num0 * a0
    t11 = num1 * a1
    t01 = num0 * num1 * c
    det = num0 * num1 * (a0 * a1 - num0 * num1 * c * c)
    return t00, t11, t01, det


def literal_selectors(spec, bl, bh, yl, yh, sl, sh, dl, dh):
    return cramer_selectors(mk(bl, bh), mk(yl, yh), mk(sl, sh), mk(dl, dh),
                            spec.get("coord") == "relative")


def _ab_coeff(n):
    """Power-series coefficients of a=(pi-phi(D))/D^2 and b=h(D)/D."""
    if n % 2 == 0:
        k = n // 2
        a = ((-1) ** k) * X.PIv / mpmath.factorial(n + 2)
        b = ((-1) ** k) * X.PIv / mpmath.factorial(n + 1)
    else:
        k = (n - 1) // 2
        a = ((-1) ** (k + 1)) * (n + 1) / mpmath.factorial(n + 2)
        b = -((-1) ** k) / mpmath.factorial(n)
    return a, b


def _series_tail(degree, radius, q=False):
    """A simple factorial majorant for the omitted analytic tail."""
    p = X.PIv
    r = X.iv.mpf([radius, radius])
    n = degree + 1
    if q:
        first = 3 * (p + n + 2) * r ** n / mpmath.factorial(n + 1)
        ratio = 2 * r / (n + 2)
    else:
        first = (p + n + 1) * r ** n / mpmath.factorial(n)
        ratio = 2 * r / (n + 1)
    return J.hi(first / (1 - ratio))


def small_gap_abq(D, degree=18):
    """Stable analytic enclosures of a, b, q=(2a-b)/D on 0<D<=3/4.

    Directly forming these removable quotients loses all useful correlation
    near D=0.  The power series is evaluated with interval arithmetic, then a
    factorial majorant encloses the entire omitted tail.
    """
    radius = abs(J.hi(D))
    if J.lo(D) <= 0 or radius > mpmath.mpf("0.75"):
        return None

    def poly(coeffs):
        out = D * 0
        for c in reversed(coeffs):
            out = out * D + c
        return out

    aa, bb = zip(*[_ab_coeff(n) for n in range(degree + 2)])
    qq = [2 * aa[n + 1] - bb[n + 1] for n in range(degree + 1)]
    ea = _series_tail(degree, radius)
    eb = _series_tail(degree, radius)
    eq = _series_tail(degree, radius, q=True)
    return (
        poly(aa[:degree + 1]) + X.iv.mpf([-ea, ea]),
        poly(bb[:degree + 1]) + X.iv.mpf([-eb, eb]),
        poly(qq) + X.iv.mpf([-eq, eq]),
    )


def _load_from_atoms(s0, s1, Z0, Z1):
    A0, _ = J.atoms(Z0)
    A1, _ = J.atoms(Z1)
    return (s0 * A0.phi + s1 * A1.phi,
            -(s0 * A0.h + s1 * A1.h),
            -(s0 * A0.sA + s1 * A1.sA))


def coordinate_loads(spec, B, Y, S, D):
    """Teacher loads at theta0 and theta1, preserving the chosen chart."""
    s0, s1 = X.masses_at(Y)
    if spec.get("coord") == "relative":
        one = _load_from_atoms(s0, s1, B + S, S)
        zero = _load_from_atoms(s0, s1, B + S + D, S + D)
    else:
        one = _load_from_atoms(s0, s1, S, S - B)
        zero = _load_from_atoms(s0, s1, S + D, S + D - B)
    return zero + one


def coordinate_rowform(spec, B, Y, S, D):
    """Division-free angle rows with affine chart correlations retained."""
    P0, A0, _, P1, A1, _ = coordinate_loads(spec, B, Y, S, D)
    AD, _ = J.atoms(D)
    M = X.PIv * X.PIv - AD.phi * AD.phi
    return ((X.PIv * P1 - AD.phi * P0) * AD.h + M * A0,
            (X.PIv * P0 - AD.phi * P1) * AD.h - M * A1)


def coordinate_mass_numerators(spec, B, Y, S, D):
    """The two exact Cramer radial numerators in the chosen chart."""
    P0, _, _, P1, _, _ = coordinate_loads(spec, B, Y, S, D)
    AD, _ = J.atoms(D)
    return X.PIv * P0 - AD.phi * P1, X.PIv * P1 - AD.phi * P0


def _gap_side(Xi, D):
    """Certified side of the gap lattice: -1 in (0,D), +1 in (D,2pi)."""
    two_pi = 2 * X.PIv
    for m in range(-4, 5):
        shifted = Xi + 2 * m * X.PIv
        if J.lo(shifted) > 0 and J.hi(shifted) < J.lo(D):
            return -1
        if J.lo(shifted) > J.hi(D) and J.hi(shifted) < J.lo(two_pi):
            return 1
    return 0


def signlaw_kills(spec, B, Y, S, D):
    """Apply the proved positive-straddle / mixed-same-side obstruction.

    This is called only after both Cramer numerators exclude zero, so an
    angle-map zero would realize a live separated critical pair.  The side
    comparisons themselves use directed interval endpoints and exact affine
    chart expressions.
    """
    yl, yh = J.lo(Y), J.hi(Y)
    if -1 < yl and yh < 0:
        positive = True
    elif 0 < yl and yh < 1:
        positive = False
    else:
        return False
    if spec.get("coord") == "relative":
        pairs = ((B + S + D, S + D), (-(B + S), -S))
    else:
        pairs = ((S + D, S + D - B), (-S, B - S))
    for left, right in pairs:
        a, b = _gap_side(left, D), _gap_side(right, D)
        if a and b and ((positive and a == b) or
                        ((not positive) and a != b)):
            return True
    return False


def collar_rowform(spec, B, Y, S, D):
    """Enclose ``(E0/D^3,E1/D^3)`` without small-gap cancellation.

    Put ``P=P(s)``, ``A=P'(s)``.  Taylor's theorem supplies independent
    enclosures ``Q=(P(s+D)-P-A*D)/D^2`` and
    ``R=(P'(s+D)-A)/D`` from the range of ``P''`` on the joining segment.
    With ``a=(pi-phi(D))/D^2``, ``b=h(D)/D`` and ``q=(2a-b)/D``, exact
    collection of the row equations gives

      E0/D^3 = b(aP-phi Q) + m R + r0 A,
      E1/D^3 = b(aP+pi Q)       + r1 A,

    where ``m=a(pi+phi)``, ``r0=pi*q+D*a*(b-a)`` and
    ``r1=-pi*q+D*a^2``.
    """
    abq = small_gap_abq(D)
    if abq is None:
        return None
    a, b, q = abq
    s0, s1 = X.masses_at(Y)
    reach = X.iv.mpf([0, J.hi(D)])
    if spec.get("coord") == "relative":
        P, A, _ = _load_from_atoms(s0, s1, B + S, S)
        _, _, P2 = _load_from_atoms(s0, s1, B + S + reach, S + reach)
    else:
        P, A, _ = _load_from_atoms(s0, s1, S, S - B)
        _, _, P2 = _load_from_atoms(s0, s1, S + reach,
                                    S + reach - B)
    Q = P2 / 2
    R = P2
    phi = X.PIv - D * D * a
    m = a * (X.PIv + phi)
    r0 = X.PIv * q + D * a * (b - a)
    r1 = -(X.PIv * q) + D * a * a
    return (b * (a * P - phi * Q) + m * R + r0 * A,
            b * (a * P + X.PIv * Q) + r1 * A)


def verdict(spec, *bx, target=None):
    """Return a valid leaf verdict, or re-derive one named ``target``.

    Search uses the default priority order.  Replay passes the stored verdict
    as ``target`` so that adding a new, earlier exclusion does not invalidate
    an older certificate whose original proof remains valid.
    """
    if target is not None and target not in VERDICTS[:-1]:
        return None
    # ``scope`` selects WHICH zeros the certificate must control.
    #
    #   "selected" (default, and what every stored artifact means): only zeros
    #       that can be local minima.  ``NEG_T00`` discharges a box by proving
    #       any zero in it is an unselected saddle, so a saddle-saddle fold is
    #       not a wall.  Consumed by the selected-count theorem.
    #
    #   "census": EVERY zero of the angle map must be fold-free, which is the
    #       ``hfold`` hypothesis of ``censusAngleMapJ_zeroCount_eq_witness``.
    #       Being unselected is then no excuse, so ``NEG_T00`` is not an
    #       admissible discharge and only a definite determinant sign or an
    #       angle-map exclusion closes a box.
    census = spec.get("scope") == "census"

    def wanted(name):
        if census and name == "NEG_T00":
            return False
        return target is None or target == name

    bl, bh, yl, yh, sl, sh, dl, dh = bx
    B, Y = mk(bl, bh), mk(yl, yh)
    S, Dm = mk(sl, sh), mk(dl, dh)
    s0, s1 = X.masses_at(Y)
    teacher = X.Teacher(B, s0, s1)
    try:
        E0, E1 = coordinate_rowform(spec, B, Y, S, Dm)
        if wanted("MAP_PLAIN") and (X.sgn(E0) or X.sgn(E1)):
            return "MAP_PLAIN"
        n0, n1 = coordinate_mass_numerators(spec, B, Y, S, Dm)
        if wanted("MAP_AUX") and X.sgn(n0) and X.sgn(n1) \
                and signlaw_kills(spec, B, Y, S, Dm):
            return "MAP_AUX"
        if wanted("MAP_CENTERED") and spec.get("coord") != "relative" and any(
                F.rowform_centered_sign(teacher, S, Dm, None)):
            return "MAP_CENTERED"
        widths = (bh - bl, yh - yl, sh - sl, dh - dl)
        if max(widths) <= spec.get("row_mv_maxw", 0.1):
            centered = rowform_mean_value(spec, *bx)
            if wanted("MAP_CENTERED") and centered is not None and (
                    X.sgn(centered[0]) or X.sgn(centered[1])):
                return "MAP_CENTERED"
        collar = collar_rowform(spec, B, Y, S, Dm)
        if wanted("MAP_COLLAR") and collar is not None and (
                X.sgn(collar[0]) or X.sgn(collar[1])):
            return "MAP_COLLAR"
        # Treat (beta,y) as interval parameters.  An empty Krawczyk image
        # excludes a zero.  More generally every zero lies in the contracted
        # angle image, so selector signs on that image classify the parent.
        if max(widths) <= spec.get("kraw_maxw", 0.1):
            kraw = coordinate_krawczyk_image(spec, *bx)
            if wanted("MAP_AUX") and kraw is not None \
                    and kraw[0] == "empty":
                return "MAP_AUX"
            if kraw is not None:
                narrow = kraw[1]
                nB, nY = B, Y
                nS, nD = mk(narrow[4], narrow[5]), mk(narrow[6], narrow[7])
                nn0, nn1 = coordinate_mass_numerators(spec, nB, nY, nS, nD)
                if wanted("MAP_AUX") and X.sgn(nn0) and X.sgn(nn1) \
                        and signlaw_kills(spec, nB, nY, nS, nD):
                    return "MAP_AUX"
                nt00, _, _, ndt = literal_selectors(spec, *narrow)
                if wanted("NEG_T00") and X.sgn(nt00) < 0:
                    return "NEG_T00"
                if wanted("NEG_DET") and X.sgn(ndt) < 0:
                    return "NEG_DET"
                if wanted("POS_DET") and X.sgn(ndt) > 0:
                    return "POS_DET"
    except Exception:
        pass

    try:
        t00, _, _, dt = literal_selectors(spec, *bx)
        if wanted("NEG_T00") and X.sgn(t00) < 0:
            return "NEG_T00"
        if wanted("NEG_DET") and X.sgn(dt) < 0:
            return "NEG_DET"
        if wanted("POS_DET") and X.sgn(dt) > 0:
            return "POS_DET"
        # AD is the expensive fallback and cannot help on coarse boxes: its
        # first-order radius is then much larger than any selector margin.
        # Delaying it is only an ordering optimization; the literal tests
        # above remain valid at every scale, and unresolved boxes are split.
        widths = (bh - bl, yh - yl, sh - sl, dh - dl)
        if max(widths) <= spec.get("mv_maxw", 0.05):
            mv = selector_mean_value(spec, *bx)
            if mv is not None:
                t00_mv, det_mv = mv
                if wanted("NEG_T00") and X.sgn(t00_mv) < 0:
                    return "NEG_T00"
                if wanted("NEG_DET") and X.sgn(det_mv) < 0:
                    return "NEG_DET"
                if wanted("POS_DET") and X.sgn(det_mv) > 0:
                    return "POS_DET"
    except Exception:
        pass
    return None


def encode(bits):
    padding = (-len(bits)) % 8
    out = bytearray()
    for i in range(0, len(bits) + padding, 8):
        byte = 0
        for k in range(8):
            byte = (byte << 1) | (bits[i + k] if i + k < len(bits) else 0)
        out.append(byte)
    return base64.b64encode(bytes(out)).decode()


def run(spec, lo_idx, hi_idx, budget):
    cells = root_cells(spec)
    bits = []
    counts = {name: 0 for name in VERDICTS}
    nodes = 0
    unresolved_volume = 0.0
    last_split = None
    started = time.time()
    for cell_idx in range(lo_idx, hi_idx):
        stack = [cells[cell_idx]]
        while stack:
            bx = stack.pop()
            nodes += 1
            if nodes > budget:
                raise RuntimeError("budget exhausted at cell %d" % cell_idx)
            if nodes % 5000 == 0:
                print("  cell %3d nodes %9d map %7d neg00 %7d negdet %7d "
                      "posdet %7d undec %d stack %6d %.0fs box=%s" %
                      (cell_idx, nodes,
                       counts["MAP_PLAIN"] + counts["MAP_CENTERED"]
                       + counts["MAP_AUX"] + counts["MAP_COLLAR"],
                       counts["NEG_T00"], counts["NEG_DET"],
                       counts["POS_DET"], counts["UNDECIDED"], len(stack),
                       time.time() - started,
                       tuple("%.6g" % float(z) for z in
                             (last_split if last_split is not None else bx))),
                      flush=True)
            result = verdict(spec, *bx)
            if result is None:
                width = max(bx[2 * i + 1] - bx[2 * i] for i in range(4))
                if width < spec["minw"]:
                    result = "UNDECIDED"
                    unresolved_volume += volume(bx)
                else:
                    last_split = bx
                    bits.append(1)
                    stack.append(child(spec, bx, 1))
                    stack.append(child(spec, bx, 0))
                    continue
            bits.append(0)
            code = CODE[result]
            bits.extend(((code >> 2) & 1, (code >> 1) & 1, code & 1))
            counts[result] += 1
    return {
        "spec": spec,
        "lo_idx": lo_idx,
        "hi_idx": hi_idx,
        "n_nodes": len(bits),
        "bits_b64": encode(bits),
        "stats": {
            "visited_nodes": nodes,
            "counts": counts,
            "unresolved_volume": unresolved_volume,
            "seconds": time.time() - started,
        },
    }


def main():
    if len(sys.argv) < 11:
        print(__doc__)
        return 2
    b0, b1, y0, y1, delta, dmax = map(float, sys.argv[1:7])
    tag = sys.argv[7]
    n_parts, lo_idx, hi_idx = map(int, sys.argv[8:11])
    budget = int(sys.argv[11]) if len(sys.argv) > 11 else 20_000_000
    minw = float(sys.argv[12]) if len(sys.argv) > 12 else 1e-8
    coord = sys.argv[13] if len(sys.argv) > 13 else "absolute"
    scope = sys.argv[14] if len(sys.argv) > 14 else "selected"
    spec = {
        "b0": b0, "b1": b1, "y0": y0, "y1": y1,
        "seam": 0.137, "delta": delta, "dmax": dmax,
        "n_parts": n_parts, "minw": minw, "mv_maxw": 0.05,
        "row_mv_maxw": 0.1, "kraw_maxw": 0.1, "coord": coord,
        "split": "relative-u-pi-endpoints-then-widest-midpoint",
    }
    # Absent from an older spec, "selected" is what it meant, so the key is
    # written only when it changes the claim.  Every stored artifact keeps
    # replaying under exactly the semantics it was produced with.
    if scope != "selected":
        spec["scope"] = scope
    if not (0 < b0 <= b1 and y0 <= y1 and 0 < delta < dmax < mpmath.pi
            and hi_idx <= n_parts and 0 <= lo_idx < hi_idx
            and coord in ("absolute", "relative")
            and scope in ("selected", "census")
            and partition_covers(spec)):
        print("SPEC REJECTED")
        return 2
    print("minimum DFS beta=[%g,%g] y=[%g,%g] D=[%g,%g] cells [%d,%d) scope=%s"
          % (b0, b1, y0, y1, delta, dmax, lo_idx, hi_idx, scope), flush=True)
    try:
        doc = run(spec, lo_idx, hi_idx, budget)
    except RuntimeError as exc:
        print(str(exc))
        return 2
    path = os.path.join(HERE, "minimum_dfs_%s_p%03d.json" % (tag, lo_idx))
    with open(path, "w") as handle:
        json.dump(doc, handle, separators=(",", ":"))
    print(json.dumps(doc["stats"], indent=1))
    print("wrote %s (%.1f KB)" %
          (os.path.basename(path), os.path.getsize(path) / 1024.0))
    if doc["stats"]["counts"]["UNDECIDED"]:
        print("INCOMPLETE: surviving boxes are PSD-wall candidates")
        return 1
    if scope == "census":
        print("COMPLETE: EVERY angle-map zero is fold-free on this region")
    else:
        print("COMPLETE: every selected zero is regular on this region")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
