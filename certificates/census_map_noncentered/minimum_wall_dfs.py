"""Three-dimensional selected-wall DFS for the noncentered census.

Teacher masses are eliminated before the sweep.  In ``(beta,s,D)``, with
``theta1=s`` and ``theta0=s+D``, a selected degenerate census zero away from
the exact fit must satisfy all three necessary conditions

  generalJAngleLocusDet = 0,
  generalJKernelTeacherSchurDetT = 0,
  generalJKernelTeacherSchurT00 >= 0.

The kernel teacher is ``(-E0(0,1), E0(1,0))``.  The last two expressions are
the same factored Cramer selectors used by ``minimum_dfs.py``, evaluated at
that teacher.  Thus a leaf is closed by a nonzero angle-locus determinant, a
nonzero kernel Schur determinant, or a negative kernel first pivot.  Only the
intersection surviving all three tests can be a genuine minimum-map wall.

The mass-free selector vanishes spuriously on the exact-fit arc; boxes wholly
inside a named tube around that arc receive ``FIT`` and are delegated to the
mass-carrying exact-fit certificate.

Usage:
  python3 minimum_wall_dfs.py b0 b1 delta dmax tag n_parts lo hi \
      [budget] [minw] [fit_radius]
"""
import base64
import json
import math
import os
import sys
import time

import mpmath

import census_cert as X
import kernel_schur as KS
import minimum_dfs as M
import noncentered as J

HERE = os.path.dirname(os.path.abspath(__file__))

VERDICTS = ("LOCUS", "REGULAR", "NEG_T00", "FIT", "UNDECIDED")
CODE = {name: i for i, name in enumerate(VERDICTS)}
NAME = {i: name for name, i in CODE.items()}


def root_box(spec):
    return (spec["b0"], spec["b1"], 0.0, float(J.hi(2 * X.PIv)),
            spec["delta"], spec["dmax"])


def root_cells(spec):
    bl, bh, sl, sh, dl, dh = root_box(spec)
    n = spec["n_parts"]
    return [(bl, bh, sl + (sh - sl) * i / n,
             sl + (sh - sl) * (i + 1) / n, dl, dh)
            for i in range(n)]


def partition_covers(spec):
    root, cells = root_box(spec), root_cells(spec)
    if not cells or cells[0][2] != root[2] or cells[-1][3] != root[3]:
        return False
    for i, cell in enumerate(cells):
        if cell[:2] != root[:2] or cell[4:] != root[4:]:
            return False
        if i and cells[i - 1][3] != cell[2]:
            return False
    return True


def split_location(bx):
    sl, sh = bx[2], bx[3]
    n0 = int(mpmath.floor(mpmath.mpf(sl) / mpmath.pi)) - 1
    n1 = int(mpmath.ceil(mpmath.mpf(sh) / mpmath.pi)) + 1
    for n in range(n0, n1 + 1):
        for cut in J.endpoints(X.iv.mpf(n) * X.PIv):
            if sl < cut < sh:
                return 1, cut
    widths = (bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4])
    axis = widths.index(max(widths))
    return axis, (bx[2 * axis] + bx[2 * axis + 1]) / 2


def child(bx, bit):
    axis, cut = split_location(bx)
    out = list(bx)
    out[2 * axis + (1 if bit == 0 else 0)] = cut
    return tuple(out)


def volume(bx):
    return ((bx[1] - bx[0]) * (bx[3] - bx[2]) *
            (bx[5] - bx[4]))


def _load(s0, s1, Z0, Z1):
    p0, h0, a0, z0 = M.atom(Z0)
    p1, h1, a1, z1 = M.atom(Z1)
    return (s0 * p0 + s1 * p1,
            -(s0 * h0 + s1 * h1),
            s0 * a0 + s1 * a1, z0 and z1)


def wall_forms(Bv, Sv, Dv):
    """Return ``(locus,T00,det,selector_smooth)`` with interval AD."""
    B = M.variable(Bv, 0)
    S = M.variable(Sv, 2)
    D = M.variable(Dv, 3)
    T1, T0 = S, S + D
    PI = M.D4(X.PIv)
    phi, h, absd, zd = M.atom(D)
    mass_det = PI * PI - phi * phi

    def rows(s0, s1):
        P0, A0, _, _ = _load(s0, s1, T0, T0 - B)
        P1, A1, _, _ = _load(s0, s1, T1, T1 - B)
        return ((PI * P1 - phi * P0) * h + mass_det * A0,
                (PI * P0 - phi * P1) * h - mass_det * A1)

    e00, e10 = rows(1, 0)
    e01, e11 = rows(0, 1)
    locus = e00 * e11 - e01 * e10
    k0, k1 = -e01, e00

    P0, _, L0, z0 = _load(k0, k1, T0, T0 - B)
    P1, _, L1, z1 = _load(k0, k1, T1, T1 - B)
    n0 = PI * P0 - phi * P1
    n1 = PI * P1 - phi * P0
    band = 2 * absd - phi
    a0 = mass_det * n1 * band - mass_det * mass_det * (2 * L0 - P0) \
        - PI * h * h * n0
    a1 = mass_det * n0 * band - mass_det * mass_det * (2 * L1 - P1) \
        - PI * h * h * n1
    c = -(mass_det * band + phi * h * h)
    t00 = n0 * a0
    det = n0 * n1 * (a0 * a1 - n0 * n1 * c * c)
    return locus, t00, det, zd and z0 and z1


def wall_forms_literal(B, S, D):
    """Cheap literal forms, used before paying for interval AD."""
    T0 = S + D
    r0 = X.sep_res_rowform(X.Teacher(B, X.thin(1), X.thin(0)),
                           T0, S, None, Dm=D)
    r1 = X.sep_res_rowform(X.Teacher(B, X.thin(0), X.thin(1)),
                           T0, S, None, Dm=D)
    locus = r0[0] * r1[1] - r1[0] * r0[1]
    det, parts = KS.kernel_teacher_schur_detT(B, T0, S)
    return locus, parts[0], det


def _collar_rows_for_masses(B, S, D, s0, s1):
    """Stable enclosures of the two angle rows divided by ``D^3``.

    This is the same exact Taylor collection as
    ``minimum_dfs.collar_rowform``, with arbitrary (possibly interval)
    teacher masses.  The student coordinates here are absolute:
    ``theta1=S`` and ``theta0=S+D``.
    """
    abq = M.small_gap_abq(D)
    if abq is None:
        return None
    a, b, q = abq
    reach = X.iv.mpf([0, J.hi(D)])
    # A high-order correlated Taylor form is decisive here.  Treating
    # Q=(P(s+D)-P-DP'(s))/D^2 and R=(P'(s+D)-P'(s))/D independently loses
    # the cancellation between their common P'' term and makes even thin
    # point intervals useless.
    deriv0 = _phi_derivatives(S, S + reach, 14)
    deriv1 = _phi_derivatives(S - B, S + reach - B, 14)
    if deriv0 is None or deriv1 is None:
        return None
    deriv = [s0 * deriv0[k] + s1 * deriv1[k]
             for k in range(len(deriv0))]
    P, A = deriv[0], deriv[1]
    phi = X.PIv - D * D * a
    mass = a * (X.PIv + phi)
    r0 = X.PIv * q + D * a * (b - a)
    r1 = -(X.PIv * q) + D * a * a
    e0 = b * a * P + r0 * A
    e1 = b * a * P + r1 * A
    for k in range(2, len(deriv)):
        power = D ** (k - 2)
        e0 += deriv[k] * power * (
            -(b * phi) / mpmath.factorial(k)
            + mass / mpmath.factorial(k - 1))
        e1 += deriv[k] * power * (
            b * X.PIv / mpmath.factorial(k))
    tail0, tail1 = _collar_tail(D, s0, s1, b, phi, mass,
                                len(deriv))
    return e0 + tail0, e1 + tail1


def _phi_derivatives(T, segment, degree):
    """Enclose ``phi^(k)(T)`` for ``0 <= k <= degree`` on one branch.

    On a half branch, ``phi=e*(sin(p)-p*cos(p))`` and ``dp/dt=-1``.
    A derivative remains ``e*((a0+a1*p)sin(p)+(b0+b1*p)cos(p))``;
    the four scalar coefficients below are advanced exactly.  Requiring the
    full Taylor segment to stay in one branch avoids making any smoothness
    assertion across a kernel kink.
    """
    # Closed half branches overlap at their endpoint because pi itself is an
    # interval.  `branch_pieces` therefore reports an extra zero-width piece
    # for [k*pi,k*pi+eps].  Select the forward branch containing the whole
    # Taylor segment instead of rejecting that legitimate one-sided Taylor
    # expansion.
    guess = int(mpmath.floor(J.mid(segment) / mpmath.pi))
    n = None
    for candidate in (guess, guess - 1, guess + 1):
        left = J.lo(X.iv.mpf(candidate) * X.PIv)
        right = J.hi(X.iv.mpf(candidate + 1) * X.PIv)
        if J.lo(segment) >= left and J.hi(segment) <= right:
            n = candidate
            break
    if n is None:
        return None
    left = J.lo(X.iv.mpf(n) * X.PIv)
    right = J.hi(X.iv.mpf(n + 1) * X.PIv)
    if J.lo(T) < left or J.hi(T) > right:
        return None
    p = J.p_of(T, n)
    e = 1 if n % 2 == 0 else -1
    sp, cp = X.iv.sin(p), X.iv.cos(p)
    a0, a1, b0, b1 = 1, 0, 0, -1
    out = []
    for _ in range(degree + 1):
        out.append(e * ((a0 + a1 * p) * sp + (b0 + b1 * p) * cp))
        # Negate the p-derivative because dp/dt=-1.
        a0, a1, b0, b1 = (-a1 + b0, b1, -a0 - b1, -a1)
    return out


def _collar_tail(D, s0, s1, b, phi, mass, first_k):
    """Symmetric rigorous majorants for the two omitted Taylor tails."""
    r = abs(J.hi(D))
    k = first_k
    # Every k-th branch derivative of phi is bounded by pi+k+1.  This follows
    # directly from the affine-sin/cos recurrence in `_phi_derivatives`.
    mass_rad = _magnitude(s0) + _magnitude(s1)
    derivative_bound = mass_rad * X.thin(J.hi(X.PIv) + k + 1)
    b_rad, phi_rad, m_rad = (_magnitude(b), _magnitude(phi),
                             _magnitude(mass))
    fact = mpmath.factorial(k)
    riv = X.iv.mpf([0, r])
    base0 = (derivative_bound * (riv ** (k - 2))
             * (b_rad * phi_rad + m_rad * k) / fact)
    base1 = (derivative_bound * (riv ** (k - 2))
             * b_rad * J.hi(X.PIv) / fact)
    ratio0 = X.iv.mpf([0, 4 * r / (k + 1)])
    ratio1 = X.iv.mpf([0, 2 * r / (k + 1)])
    err0 = J.hi(base0 / (1 - ratio0))
    err1 = J.hi(base1 / (1 - ratio1))
    return X.iv.mpf([-err0, err0]), X.iv.mpf([-err1, err1])


def _magnitude(z):
    """Outward interval upper bound for ``abs(z)`` (never a float cast)."""
    return X.iv.mpf([0, max(abs(J.lo(z)), abs(J.hi(z)))])


def wall_forms_normalized(B, S, D):
    """Small-gap wall forms with all forced powers of ``D`` removed.

    Let ``Eij`` be the two angle rows for the two unit teachers.  Each row is
    ``D^3 * eij``.  We use the normalized angle-kernel teacher

        ``k = (-e01, e00)``.

    The unnormalized kernel is ``D^3*k``.  Homogeneity therefore gives the
    same selector signs.  For this normalized teacher, write

        ``N_i = D*n_i``, ``A_i = D^3*a_i``, ``C = D^2*c``.

    Then the returned quantities are exactly the angle-locus determinant
    divided by ``D^6``, the first kernel selector divided by ``D^4`` and its
    determinant divided by ``D^8`` for the normalized kernel teacher.  The
    corresponding forms for the original angle kernel differ by the positive
    factors ``D^6`` and ``D^12`` from mass homogeneity.  Thus all three signs
    and zero sets are unchanged for ``D>0``.

    No interval division is used.  Every removable quotient is evaluated by
    the power-series enclosure in ``small_gap_abq`` and Taylor's theorem.
    """
    rows0 = _collar_rows_for_masses(B, S, D, D * 0 + 1, D * 0)
    rows1 = _collar_rows_for_masses(B, S, D, D * 0, D * 0 + 1)
    if rows0 is None or rows1 is None:
        return None
    e00, e10 = rows0
    e01, e11 = rows1
    locus = e00 * e11 - e01 * e10
    primary = _selector_for_normalized_kernel(B, S, D, -e01, e00)
    alternate = _selector_for_normalized_kernel(B, S, D, -e11, e10)
    if primary is None or alternate is None:
        return None
    # The Boolean flags certify that the row from which the kernel was built
    # cannot be the zero row anywhere in the box.  A selector verdict is used
    # only with its corresponding flag.  The alternate row removes the
    # spurious zero kernel on torque-dead faces such as theta1=0.
    first_alive = bool(X.sgn(e00) or X.sgn(e01))
    second_alive = bool(X.sgn(e10) or X.sgn(e11))
    return (locus, primary[0], primary[1], alternate[0], alternate[1],
            first_alive, second_alive)


def _selector_for_normalized_kernel(B, S, D, k0, k1):
    """Normalized Cramer ``(T00,det)`` for one angle-row kernel."""
    reach = X.iv.mpf([0, J.hi(D)])
    d0 = _phi_derivatives(S, S + reach, 14)
    d1 = _phi_derivatives(S - B, S + reach - B, 14)
    if d0 is None or d1 is None:
        return None
    deriv = [k0 * d0[j] + k1 * d1[j] for j in range(len(d0))]
    P, A = deriv[0], deriv[1]
    abq = M.small_gap_abq(D)
    a, b, _ = abq
    phi = X.PIv - D * D * a
    mass = a * (X.PIv + phi)       # (pi^2-phi^2)/D^2
    # Correlated Taylor series for N0/D and N1/D.  Independent remainder
    # quotients lose the leading cancellation in the Schur pivot.
    n0 = X.PIv * A + D * a * P
    n1 = -X.PIv * A + D * a * P + D * D * a * A
    for j in range(2, len(deriv)):
        fact = mpmath.factorial(j)
        n0 += X.PIv * deriv[j] * D ** (j - 1) / fact
        n1 += deriv[j] * (
            -(X.PIv * D ** (j - 1)) + a * D ** (j + 1)) / fact
    ntail0, ntail1 = _numerator_tails(D, k0, k1, a, len(deriv))
    n0 += ntail0
    n1 += ntail1
    P1 = P
    atom1, _ = J.atoms(S)
    atom1b, _ = J.atoms(S - B)
    atom0, _ = J.atoms(S + D)
    atom0b, _ = J.atoms(S + D - B)
    P0 = k0 * atom0.phi + k1 * atom0b.phi

    L1 = k0 * atom1.asin + k1 * atom1b.asin
    L0 = k0 * atom0.asin + k1 * atom0b.asin
    atomd, _ = J.atoms(D)
    band = 2 * atomd.asin - phi
    c = -(mass * band + phi * b * b)
    a0 = (mass * n1 * band - D * mass * mass * (2 * L0 - P0)
          - X.PIv * b * b * n0)
    a1 = (mass * n0 * band - D * mass * mass * (2 * L1 - P1)
          - X.PIv * b * b * n1)
    t00 = n0 * a0
    det = n0 * n1 * (a0 * a1 - n0 * n1 * c * c)
    return t00, det


def _numerator_tails(D, k0, k1, a, first_k):
    """Majorants for the omitted normalized Cramer-numerator series."""
    r = abs(J.hi(D))
    k = first_k
    mass_rad = _magnitude(k0) + _magnitude(k1)
    derivative_bound = mass_rad * X.thin(J.hi(X.PIv) + k + 1)
    fact = mpmath.factorial(k)
    common = derivative_bound * X.iv.mpf([0, r]) ** (k - 1) / fact
    base0 = common * J.hi(X.PIv)
    base1 = common * (X.thin(J.hi(X.PIv)) + _magnitude(a) * r * r)
    ratio = X.iv.mpf([0, 2 * r / (k + 1)])
    err0 = J.hi(base0 / (1 - ratio))
    err1 = J.hi(base1 / (1 - ratio))
    return X.iv.mpf([-err0, err0]), X.iv.mpf([-err1, err1])


def mean_values(boxes, full):
    mids = tuple(X.thin(J.mid(z)) for z in boxes)
    center = wall_forms(*mids)

    def form(which):
        q = center[which].v
        for axis in (0, 2, 3):
            source = 0 if axis == 0 else axis - 1
            q += full[which].g[axis] * (boxes[source] - mids[source])
        return q

    return form(0), form(1), form(2)


def inside_fit(spec, bx):
    bl, bh, sl, sh, dl, dh = bx
    r = spec["fit_radius"]
    two_pi = float(J.hi(2 * X.PIv))
    seam_side = (0 <= sl and sh <= r) or (two_pi - r <= sl and sh <= two_pi)
    return seam_side and dl - bh >= -r and dh - bl <= r


def verdict(spec, *bx):
    if inside_fit(spec, bx):
        return "FIT"
    boxes = (M.mk(bx[0], bx[1]), M.mk(bx[2], bx[3]), M.mk(bx[4], bx[5]))
    try:
        normalized = wall_forms_normalized(*boxes)
        if normalized is not None:
            if X.sgn(normalized[0]):
                return "LOCUS"
            for t00, det, alive in (
                    (normalized[1], normalized[2], normalized[5]),
                    (normalized[3], normalized[4], normalized[6])):
                if not alive:
                    continue
                if X.sgn(det):
                    return "REGULAR"
                if X.sgn(t00) < 0:
                    return "NEG_T00"
        literal = wall_forms_literal(*boxes)
        locus = literal[0]
        if X.sgn(locus):
            return "LOCUS"
        if X.sgn(literal[2]):
            return "REGULAR"
        if X.sgn(literal[1]) < 0:
            return "NEG_T00"
        widths = (bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4])
        if max(widths) <= spec.get("mv_maxw", 0.15):
            full = wall_forms(*boxes)
            locus_mv, t00_mv, det_mv = mean_values(boxes, full)
            if X.sgn(locus_mv):
                return "LOCUS"
            if full[3]:
                if X.sgn(det_mv):
                    return "REGULAR"
                if X.sgn(t00_mv) < 0:
                    return "NEG_T00"
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
    candidates = []
    nodes = 0
    started = time.time()
    for idx in range(lo_idx, hi_idx):
        stack = [cells[idx]]
        while stack:
            bx = stack.pop()
            nodes += 1
            if nodes > budget:
                raise RuntimeError("budget exhausted at cell %d" % idx)
            if nodes % 5000 == 0:
                print("cell %d nodes %d counts=%s stack=%d %.0fs box=%s" %
                      (idx, nodes, counts, len(stack), time.time() - started,
                       tuple("%.6g" % float(z) for z in bx)), flush=True)
            result = verdict(spec, *bx)
            if result is None:
                width = max(bx[1] - bx[0], bx[3] - bx[2], bx[5] - bx[4])
                if width < spec["minw"]:
                    result = "UNDECIDED"
                    if len(candidates) < 1000:
                        candidates.append(tuple(float(z) for z in bx))
                else:
                    bits.append(1)
                    stack.append(child(bx, 1))
                    stack.append(child(bx, 0))
                    continue
            bits.append(0)
            code = CODE[result]
            bits.extend(((code >> 2) & 1, (code >> 1) & 1, code & 1))
            counts[result] += 1
    return {"spec": spec, "lo_idx": lo_idx, "hi_idx": hi_idx,
            "n_nodes": len(bits), "bits_b64": encode(bits),
            "stats": {"visited_nodes": nodes, "counts": counts,
                      "seconds": time.time() - started,
                      "candidate_boxes": candidates}}


def main():
    if len(sys.argv) < 9:
        print(__doc__)
        return 2
    b0, b1, delta, dmax = map(float, sys.argv[1:5])
    tag = sys.argv[5]
    n_parts, lo_idx, hi_idx = map(int, sys.argv[6:9])
    budget = int(sys.argv[9]) if len(sys.argv) > 9 else 5_000_000
    minw = float(sys.argv[10]) if len(sys.argv) > 10 else 1e-5
    fit_radius = float(sys.argv[11]) if len(sys.argv) > 11 else 0.01
    spec = {"b0": b0, "b1": b1, "delta": delta, "dmax": dmax,
            "n_parts": n_parts, "minw": minw, "fit_radius": fit_radius,
            "mv_maxw": 0.15, "split": "s-pi-endpoints-then-widest"}
    if not (0 < b0 <= b1 < math.pi and 0 < delta < dmax < math.pi
            and 0 <= lo_idx < hi_idx <= n_parts and partition_covers(spec)):
        print("SPEC REJECTED")
        return 2
    print("minimum wall DFS beta=[%g,%g] cells=[%d,%d)" %
          (b0, b1, lo_idx, hi_idx), flush=True)
    try:
        doc = run(spec, lo_idx, hi_idx, budget)
    except RuntimeError as exc:
        print(exc)
        return 2
    path = os.path.join(HERE, "minimum_wall_%s_p%03d.json" % (tag, lo_idx))
    with open(path, "w") as handle:
        json.dump(doc, handle, separators=(",", ":"))
    print(json.dumps(doc["stats"], indent=1))
    print("wrote", os.path.basename(path))
    return 1 if doc["stats"]["counts"]["UNDECIDED"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
