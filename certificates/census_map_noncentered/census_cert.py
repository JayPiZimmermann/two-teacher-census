"""
The CERTIFIED CENSUS of a single noncentered teacher (beta, s0, s1).

Two halves, both certified in mpmath.iv:

COINCIDENT half.  The torque A(t) = -(s0 h(t) + s1 h(t-beta)) is C^1 (h' =
slopeAtomJ is continuous), so its roots are isolated by interval Newton with a
certified-complete branch and bound over t in [0, 2pi].  At each root the mass
direction is (s0,s1) = kappa*(-h(t-beta), h(t)) for some real kappa != 0, hence
      P   = -kappa*Wpot ,  W = -kappa*Wwgt ,  tau = P - 2W = +kappa*Wtau
the last by the identity Wtau + Wpot = 2 Wwgt (validate.py (d), sympy-exact for
this kernel).  kappa^2 cancels in widgets.js coincidenceTypeAt, so
      a trap is present at the root   <=>   Wtau*Wpot < 0
      its label is @positive          <=>   Wtau*Wwgt > 0     (else @mixed)
exactly as in the centered arrangement.

SEPARATED half.  Two students at th0 != th1 with weights c0, c1.  The RADIAL
rows solve the weights in the WELL-CONDITIONED direction,
      massDet = kap^2 - phi(D)^2 ,  D = th0 - th1,  kap = pi
      c0 = (kap P(th0) - phi(D) P(th1)) / massDet
      c1 = (kap P(th1) - phi(D) P(th0)) / massDet
and the two ANGULAR rows are the residual
      G0 = c1 h(D) + A(th0) ,   G1 = c0 h(D) - A(th1).
massDet = pi^2 - phiJ(D)^2 vanishes ONLY on D = 0 (mod 2pi) -- phiJ has range
[0, pi] with phiJ = pi only at D = 0 -- so this direction is singular only at
student coincidence, where the coincident rows apply instead.  (widgets.js
separatedScan divides by h(D) instead, which is singular on D = 0 AND D = pi
and is guarded by |h(D)| > 5e-3.)
The zeros of (G0,G1) are isolated and are enumerated by a certified branch and
bound with a 2x2 KRAWCZYK test for existence and uniqueness.
Each solution's type is widgets.js schurLabel, evaluated in interval arithmetic.
"""
import math

import mpmath
from mpmath import iv

import noncentered as J

PREC = 160
J.set_prec(PREC)
PIv = J.PI_IV()


def thin(x):
    return iv.mpf([x, x])


def inter(x, y):
    a1, b1 = J.endpoints(x)
    a2, b2 = J.endpoints(y)
    a, b = max(a1, a2), min(b1, b2)
    if a > b:
        return None
    return iv.mpf([a, b])


def sgn(x):
    a, b = J.endpoints(x)
    if a > 0:
        return 1
    if b < 0:
        return -1
    return 0


# --------------------------------------------------------------------------
# teacher loads
# --------------------------------------------------------------------------

class Teacher(object):
    def __init__(self, beta, s0, s1):
        self.B = beta
        self.s0 = s0
        self.s1 = s1

    def load(self, T, branches=None):
        """(P, A, Ad) at angle T: P = s0 phi(t)+s1 phi(t-beta), A = P' ,
        Ad = A'."""
        if branches is None:
            A0, _ = J.atoms(T)
            A1, _ = J.atoms(T - self.B)
        else:
            A0 = J.atoms_on_branch(T, branches[0])
            A1 = J.atoms_on_branch(T - self.B, branches[1])
        P = self.s0 * A0.phi + self.s1 * A1.phi
        A = -(self.s0 * A0.h + self.s1 * A1.h)
        Ad = -(self.s0 * A0.sA + self.s1 * A1.sA)
        return P, A, Ad


def masses_at(y):
    """widgets.js massesAt: psi = (y+1)pi/2, (s0,s1) = (sin psi, cos psi)."""
    psi = (y + 1) * PIv / 2
    return iv.sin(psi), iv.cos(psi)


# --------------------------------------------------------------------------
# COINCIDENT half: certified-complete torque roots + sign chart
# --------------------------------------------------------------------------

SPLIT = mpmath.mpf("0.5173648177666931")   # off-centre split (see ARRANGEMENT)


def torque_roots(Tc, tol=mpmath.mpf("1e-28"), maxdepth=90):
    """Certified-complete isolation of the roots of A(t) on [0, 2pi].

    Returns (roots, undecided).  Each root is a dict with the enclosure and the
    certified signs of Wtau, Wpot, Wwgt."""
    two_pi = J.hi(2 * PIv)
    undecided = []
    found = []
    stack = [(mpmath.mpf(0), two_pi, 0)]
    while stack:
        a, b, dep = stack.pop()
        T = iv.mpf([a, b])
        P, A, Ad = Tc.load(T)
        if sgn(A) != 0:
            continue
        if b - a <= tol:
            found.append((a, b))
            continue
        if sgn(Ad) != 0:
            sa = sgn(Tc.load(thin(a))[1])
            sb = sgn(Tc.load(thin(b))[1])
            if sa != 0 and sb != 0:
                if sa == sb:
                    continue
                # exactly one root: bisect then interval Newton
                lo_, hi_ = a, b
                for _ in range(80):
                    m = (lo_ + hi_) / 2
                    sm = sgn(Tc.load(thin(m))[1])
                    if sm == 0:
                        break
                    if sm == sa:
                        lo_ = m
                    else:
                        hi_ = m
                    if hi_ - lo_ <= tol:
                        break
                Tr = iv.mpf([lo_, hi_])
                for _ in range(20):
                    m = thin(J.mid(Tr))
                    dd = Tc.load(Tr)[2]
                    if sgn(dd) == 0:
                        break
                    N = m - Tc.load(m)[1] / dd
                    NN = inter(N, Tr)
                    if NN is None or J.width(NN) >= J.width(Tr):
                        break
                    Tr = NN
                found.append(J.endpoints(Tr))
                continue
        if dep >= maxdepth:
            undecided.append((float(a), float(b)))
            continue
        m = a + (b - a) * SPLIT
        stack.append((a, m, dep + 1))
        stack.append((m, b, dep + 1))
    found.sort()
    # merge touching enclosures (a root exactly on a subdivision point)
    merged = []
    for a, b in found:
        if merged and a <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], b))
        else:
            merged.append((a, b))
    # t = 0 and t = 2pi are the same point
    if len(merged) >= 2 and merged[0][0] <= 0 and merged[-1][1] >= two_pi:
        merged = merged[:-1]
    out = []
    for a, b in merged:
        T = iv.mpf([a, b])
        Dts = J.dets_iv(Tc.B, T)
        rec = {
            "t": [mpmath.nstr(a, 25), mpmath.nstr(b, 25)],
            "t_mid": float((a + b) / 2),
            "sgn_Wtau": sgn(Dts["Wtau"][0]),
            "sgn_Wpot": sgn(Dts["Wpot"][0]),
            "sgn_Wwgt": sgn(Dts["Wwgt"][0]),
        }
        st, sp, sw = rec["sgn_Wtau"], rec["sgn_Wpot"], rec["sgn_Wwgt"]
        if st == 0 or sp == 0 or sw == 0:
            rec["label"] = None
        elif st * sp < 0:
            rec["label"] = ("coincident:trap@positive" if st * sw > 0
                            else "coincident:trap@mixed")
        else:
            rec["label"] = ""          # saddle, contributes nothing
        out.append(rec)
    return out, undecided


# --------------------------------------------------------------------------
# SEPARATED half: certified-complete 2-D enumeration
# --------------------------------------------------------------------------

def sep_res(Tc, T0, T1, br=None, want_jac=False, Dm=None):
    """(G0, G1, c0, c1, massDet) and optionally the 2x2 Jacobian.

    Dm is the D = th0 - th1 box.  It MUST be passed when the caller knows it
    exactly: computing T0 - T1 inside interval arithmetic widens D by twice the
    width of T1 (T0 is built as T1 + D), which is a pure dependency artefact and
    was large enough to make massDet straddle zero on perfectly good boxes."""
    if Dm is None:
        Dm = T0 - T1
    if br is None:
        AD, _ = J.atoms(Dm)
        b0 = b1 = None
    else:
        AD = J.atoms_on_branch(Dm, br["nD"])
        b0 = (br["n0"], br["n0b"])
        b1 = (br["n1"], br["n1b"])
    hD, phiD, sAD = AD.h, AD.phi, AD.sA
    kap = PIv
    M = kap * kap - phiD * phiD
    if J.contains_zero(M):
        return None
    P0, A0, Ad0 = Tc.load(T0, b0)
    P1, A1, Ad1 = Tc.load(T1, b1)
    c0 = (kap * P0 - phiD * P1) / M
    c1 = (kap * P1 - phiD * P0) / M
    G0 = c1 * hD + A0
    G1 = c0 * hD - A1
    if not want_jac:
        return G0, G1, c0, c1, M, hD, phiD, sAD
    # D = th0 - th1 and phi' = -h, so d(phiD)/d(th0) = -hD and
    # d(phiD)/d(th1) = +hD, while dM/d(th0) = q = 2 phiD hD = -dM/d(th1).
    # SIGN FIX 2026-08-08: dc1_1 carried `+ hD * P0` where the chain rule
    # gives `- hD * P0`.  That made the INTERVAL j01 an invalid enclosure of
    # dG0/d(th1) (measured: it missed the true derivative range in 188 of 300
    # random boxes) and hence made the Krawczyk test return spurious
    # `unique` verdicts -- certifying separated families that do not exist.
    # Regression guard: sep_jacobian_audit.py.
    q = 2 * phiD * hD                      # dM/dD
    dc0_0 = (kap * A0 + hD * P1) / M - c0 * q / M
    dc0_1 = (-(hD * P1) - phiD * A1) / M + c0 * q / M
    dc1_0 = (hD * P0 - phiD * A0) / M - c1 * q / M
    dc1_1 = (kap * A1 - hD * P0) / M + c1 * q / M
    j00 = dc1_0 * hD + c1 * sAD + Ad0
    j01 = dc1_1 * hD - c1 * sAD
    j10 = dc0_0 * hD + c0 * sAD
    j11 = dc0_1 * hD - c0 * sAD - Ad1
    return G0, G1, c0, c1, M, hD, phiD, sAD, (j00, j01, j10, j11)


def sep_res_rowform(Tc, T0, T1, br=None, Dm=None):
    """(E0, E1), the DIVISION-FREE row form of the two eliminated angle
    equations.  `E_i = massDet(D) * G_i` with the `G_i` of `sep_res`, so the
    zero sets coincide and a sign-definite `E_i` excludes the box exactly as a
    sign-definite `G_i` does.

    WHY IT IS A DIFFERENT TEST.  `sep_res` forms the Cramer weights by DIVIDING
    by `M = pi^2 - phi(D)^2` and then multiplies by `hD`; in interval
    arithmetic the division inflates by the relative width of `M` and the later
    multiplication cannot take it back.  The row form never divides, so its
    enclosure of the same quantity is at least as tight and is usually much
    tighter on boxes where `D` is wide.  It also survives `M` straddling zero,
    where `sep_res` must give up and return `None`.

    THE IDENTITY IS THE TREE'S, not a rederivation.  With
    `num(D,x) = (pi*phi(x-D) - phi(D)*phi(x))*h(D) - (pi^2-phi(D)^2)*h(x)`
    (`SeparatedCount/ScalarColumnJ.separatedNumJ`), the landed
    `generalJAngleEq0_linear` + `generalJAngleEq0_unitTeacher0/1` give
    `E0 = s0*num(D,th0) + s1*num(D,th0-beta)`, and collecting the teacher
    masses turns that into the two lines below; `generalJAngleEq1_*` do the
    same for `E1`.  `validate_rowform.py` is the numerical gate.
    """
    if Dm is None:
        Dm = T0 - T1
    if br is None:
        AD, _ = J.atoms(Dm)
        b0 = b1 = None
    else:
        AD = J.atoms_on_branch(Dm, br["nD"])
        b0 = (br["n0"], br["n0b"])
        b1 = (br["n1"], br["n1b"])
    hD, phiD = AD.h, AD.phi
    kap = PIv
    M = kap * kap - phiD * phiD
    P0, A0, _ = Tc.load(T0, b0)
    P1, A1, _ = Tc.load(T1, b1)
    E0 = (kap * P1 - phiD * P0) * hD + M * A0
    E1 = (kap * P0 - phiD * P1) * hD - M * A1
    return E0, E1


def sep_res_rowform_jac(Tc, T0, T1, br=None, Dm=None):
    """`(E0, E1, dE0/ds, dE0/dD, dE1/ds, dE1/dD)` -- the row form AND its
    derivative, both division-free, in the `(s, D)` coordinates
    `th0 = s + D`, `th1 = s`.

    WHY THIS EXISTS.  Measured at the F1 witness on 200 plateau-scale boxes
    that survive every current test: the CRAMER Jacobian entry `dG0/dD` is
    enclosed **5041x wider than its true range** (median; p90 2.9e5), against
    512x for the residual itself.  Each entry divides by `M` three times and
    then subtracts two same-sized quotients (`(...)/M - c*q/M`), which is a
    textbook dependency amplifier.  Everything that consumes the Jacobian
    inherits it: the mean-value form loses its second-order advantage, the
    monotonicity test almost never sees a sign-definite partial (6 of 308
    boxes), and Krawczyk contracts far less than it should.

    Differentiating `E = M*G` instead never divides at all.  With
    `phi' = -h`, `h' = sA`, `P' = A`, `A' = Ad`, and `dM/dD = 2*phi(D)*h(D)`:

        dE0/ds = (kap*A1 - phiD*A0)*hD + M*Ad0
        dE0/dD = (hD*P0 - phiD*A0)*hD + (kap*P1 - phiD*P0)*sAD
                 + 2*phiD*hD*A0 + M*Ad0
        dE1/ds = (kap*A0 - phiD*A1)*hD - M*Ad1
        dE1/dD = (kap*A0 + hD*P1)*hD + (kap*P0 - phiD*P1)*sAD
                 - 2*phiD*hD*A1

    `validate_rowform.py` gates these two ways: against finite differences,
    and against the audited Cramer Jacobian through the exact relations
    `dE/ds = M*dG/ds` and `dE0/dD = 2*phiD*hD*G0 + M*dG0/dD`.  A derivative
    enclosure that is merely INVALID is the failure this campaign has already
    paid for once (the 2026-08-08 sign fix), so it is gated, not trusted.
    """
    if Dm is None:
        Dm = T0 - T1
    if br is None:
        AD, _ = J.atoms(Dm)
        b0 = b1 = None
    else:
        AD = J.atoms_on_branch(Dm, br["nD"])
        b0 = (br["n0"], br["n0b"])
        b1 = (br["n1"], br["n1b"])
    hD, phiD, sAD = AD.h, AD.phi, AD.sA
    kap = PIv
    M = kap * kap - phiD * phiD
    P0, A0, Ad0 = Tc.load(T0, b0)
    P1, A1, Ad1 = Tc.load(T1, b1)
    U = kap * P1 - phiD * P0
    V = kap * P0 - phiD * P1
    E0 = U * hD + M * A0
    E1 = V * hD - M * A1
    q = 2 * phiD * hD
    dE0_s = (kap * A1 - phiD * A0) * hD + M * Ad0
    dE0_D = (hD * P0 - phiD * A0) * hD + U * sAD + q * A0 + M * Ad0
    dE1_s = (kap * A0 - phiD * A1) * hD - M * Ad1
    dE1_D = (kap * A0 + hD * P1) * hD + V * sAD - q * A1
    return E0, E1, dE0_s, dE0_D, dE1_s, dE1_D


def krawczyk(Tc, T0, T1, br, Dm=None):
    """2x2 Krawczyk test. Returns 'unique', 'empty' or None."""
    r = sep_res(Tc, T0, T1, br, want_jac=True, Dm=Dm)
    if r is None:
        return None
    j00, j01, j10, j11 = r[8]
    x0, x1 = thin(J.mid(T0)), thin(J.mid(T1))
    rc = sep_res(Tc, x0, x1, br, want_jac=True)
    if rc is None:
        return None
    m00, m01, m10, m11 = rc[8]
    det = J.mid(m00) * J.mid(m11) - J.mid(m01) * J.mid(m10)
    if det == 0:
        return None
    y00 = J.mid(m11) / det
    y01 = -J.mid(m01) / det
    y10 = -J.mid(m10) / det
    y11 = J.mid(m00) / det
    Y = (thin(y00), thin(y01), thin(y10), thin(y11))
    F0, F1 = rc[0], rc[1]
    # I - Y*J
    a00 = 1 - (Y[0] * j00 + Y[1] * j10)
    a01 = -(Y[0] * j01 + Y[1] * j11)
    a10 = -(Y[2] * j00 + Y[3] * j10)
    a11 = 1 - (Y[2] * j01 + Y[3] * j11)
    r0 = T0 - x0
    r1 = T1 - x1
    K0 = x0 - (Y[0] * F0 + Y[1] * F1) + (a00 * r0 + a01 * r1)
    K1 = x1 - (Y[2] * F0 + Y[3] * F1) + (a10 * r0 + a11 * r1)
    if inter(K0, T0) is None or inter(K1, T1) is None:
        return "empty"
    a, b = J.endpoints(T0)
    ka, kb = J.endpoints(K0)
    c, d = J.endpoints(T1)
    la, lb = J.endpoints(K1)
    if a < ka and kb < b and c < la and lb < d:
        return "unique"
    return None


def sep_refine(Tc, T0, T1, br, iters=40, Dm=None):
    for _ in range(iters):
        r = sep_res(Tc, T0, T1, br, want_jac=True, Dm=Dm)
        if r is None:
            break
        j00, j01, j10, j11 = r[8]
        x0, x1 = thin(J.mid(T0)), thin(J.mid(T1))
        rc = sep_res(Tc, x0, x1, br, want_jac=True)
        if rc is None:
            break
        m00, m01, m10, m11 = rc[8]
        det = J.mid(m00) * J.mid(m11) - J.mid(m01) * J.mid(m10)
        if det == 0:
            break
        Y = (thin(J.mid(m11) / det), thin(-J.mid(m01) / det),
             thin(-J.mid(m10) / det), thin(J.mid(m00) / det))
        F0, F1 = rc[0], rc[1]
        a00 = 1 - (Y[0] * j00 + Y[1] * j10)
        a01 = -(Y[0] * j01 + Y[1] * j11)
        a10 = -(Y[2] * j00 + Y[3] * j10)
        a11 = 1 - (Y[2] * j01 + Y[3] * j11)
        K0 = x0 - (Y[0] * F0 + Y[1] * F1) + (a00 * (T0 - x0) + a01 * (T1 - x1))
        K1 = x1 - (Y[2] * F0 + Y[3] * F1) + (a10 * (T0 - x0) + a11 * (T1 - x1))
        N0, N1 = inter(K0, T0), inter(K1, T1)
        if N0 is None or N1 is None:
            break
        if J.width(N0) >= J.width(T0) and J.width(N1) >= J.width(T1):
            break
        T0, T1 = N0, N1
        Dm = None
    return T0, T1


def schur_label(Tc, T0, T1, c0, c1):
    """widgets.js schurLabel in interval arithmetic; returns
    ('spurious'|'saddle'|'undecided', signs)."""
    Dm = T0 - T1
    AD, _ = J.atoms(Dm)
    hD, phiD, sAD = AD.h, AD.phi, AD.sA
    kap = PIv
    P0, A0, _ = Tc.load(T0)
    P1, A1, _ = Tc.load(T1)
    # Tau(K,T,t) = Pot - 2 Wgt ;  Wgt = s0|sin t| + s1|sin(t-beta)|
    At0, _ = J.atoms(T0)
    At0b, _ = J.atoms(T0 - Tc.B)
    At1, _ = J.atoms(T1)
    At1b, _ = J.atoms(T1 - Tc.B)
    W0 = Tc.s0 * At0.asin + Tc.s1 * At0b.asin
    W1 = Tc.s0 * At1.asin + Tc.s1 * At1b.asin
    tau0 = P0 - 2 * W0
    tau1 = P1 - 2 * W1
    Hp = sAD                                     # K.dH(D)
    A00 = -(c0 * c1 * Hp) + c0 * tau0
    A11 = -(c0 * c1 * Hp) + c1 * tau1
    A01 = c0 * c1 * Hp
    M = kap * kap - phiD * phiD
    T00 = M * A00 - kap * hD * hD * c0 * c0
    T11 = M * A11 - kap * hD * hD * c1 * c1
    T01 = M * A01 - phiD * hD * hD * c0 * c1
    detT = T00 * T11 - T01 * T01
    s00, s11, sdet = sgn(T00), sgn(T11), sgn(detT)
    if s00 > 0 and sdet > 0:
        lab = "spurious"
    elif s00 < 0 or s11 < 0 or sdet < 0:
        lab = "saddle"
    elif s00 == 0 or sdet == 0:
        lab = "undecided_interval"
    else:
        lab = "undecided_second_order"
    return lab, (s00, s11, sdet)


def _krawczyk_sD(Tc, T1, Dm, br):
    """Krawczyk in the (s, D) = (th1, th0-th1) coordinates, on a genuine box.

    th0 = s + D is exact as a range, so the cell IS a box in (s,D) and a
    "unique" verdict is a statement about that cell alone -- two disjoint cells
    can never claim the same solution.  Chain rule:
        d/ds = d/dth0 + d/dth1 ,   d/dD = d/dth0 ."""
    T0 = T1 + Dm
    r = sep_res(Tc, T0, T1, br, want_jac=True, Dm=Dm)
    if r is None:
        return None, None, None
    j00, j01, j10, j11 = r[8]
    A00, A01 = j00 + j01, j00
    A10, A11 = j10 + j11, j10
    sm, dm = thin(J.mid(T1)), thin(J.mid(Dm))
    rc = sep_res(Tc, sm + dm, sm, br, want_jac=True, Dm=dm)
    if rc is None:
        return None, None, None
    m00, m01, m10, m11 = rc[8]
    c00, c01 = J.mid(m00 + m01), J.mid(m00)
    c10, c11 = J.mid(m10 + m11), J.mid(m10)
    det = c00 * c11 - c01 * c10
    if det == 0:
        return None, None, None
    Y = (thin(c11 / det), thin(-c01 / det), thin(-c10 / det), thin(c00 / det))
    F0, F1 = rc[0], rc[1]
    b00 = 1 - (Y[0] * A00 + Y[1] * A10)
    b01 = -(Y[0] * A01 + Y[1] * A11)
    b10 = -(Y[2] * A00 + Y[3] * A10)
    b11 = 1 - (Y[2] * A01 + Y[3] * A11)
    rs, rd = T1 - sm, Dm - dm
    K0 = sm - (Y[0] * F0 + Y[1] * F1) + (b00 * rs + b01 * rd)
    K1 = dm - (Y[2] * F0 + Y[3] * F1) + (b10 * rs + b11 * rd)
    N0, N1 = inter(K0, T1), inter(K1, Dm)
    if N0 is None or N1 is None:
        return "empty", None, None
    a, b = J.endpoints(T1)
    ka, kb = J.endpoints(K0)
    c, d = J.endpoints(Dm)
    la, lb = J.endpoints(K1)
    if a < ka and kb < b and c < la and lb < d:
        return "unique", N0, N1
    return None, None, None


def separated_families(Tc, delta=mpmath.mpf("0.05"), minw=mpmath.mpf("2e-5"),
                       budget=120000, cell=0.35):
    """Certified-complete enumeration of the separated families of ONE teacher.

    Domain: s = th1 in [0,2pi], D = th0-th1 in [delta, pi].  Because the two
    students are UNLABELLED, every unordered pair {th0,th1} with th0 != th1
    appears in that domain EXACTLY ONCE (of D and 2pi-D exactly one is <= pi),
    so no swap de-duplication is needed and none is done.
    The strip D in [0, delta) is the COINCIDENCE COLLAR: massDet = pi^2-phiJ(D)^2
    vanishes only at D = 0 (mod 2pi), the weights do not solve there and the
    coincident rows apply instead.  Its area is reported.
    D = pi (the ANTIPODAL stratum, h(D) = 0) is INSIDE the certified domain:
    the residual used here never divides by h(D), unlike widgets.js
    separatedScan which guards |h(D)| > 5e-3 and cannot see it."""
    two_pi = J.hi(2 * PIv)
    dpi = J.hi(PIv)
    cuts = {mpmath.mpf(0), two_pi, dpi}
    a, b = J.endpoints(Tc.B)
    for z in (a, b, a + dpi, b + dpi, a - dpi, b - dpi):
        if 0 < z < two_pi:
            cuts.add(z)
    cuts = sorted(cuts)
    sf = []
    for i in range(len(cuts) - 1):
        u, v = cuts[i], cuts[i + 1]
        n = max(1, int(math.ceil(float(v - u) / cell)))
        for k in range(n):
            sf.append((u + (v - u) * k / n, u + (v - u) * (k + 1) / n))
    dcuts = sorted({delta, dpi} | {z for z in (a, b, a + dpi, b + dpi,
                                               two_pi - b, two_pi - a)
                                   if delta < z < dpi})
    df = []
    for i in range(len(dcuts) - 1):
        u, v = dcuts[i], dcuts[i + 1]
        n = max(1, int(math.ceil(float(v - u) / cell)))
        for k in range(n):
            df.append((u + (v - u) * k / n, u + (v - u) * (k + 1) / n))

    stack = [(x[0], x[1], z[0], z[1]) for x in sf for z in df]
    sols = []
    undecided = []
    steps = 0
    while stack:
        steps += 1
        if steps > budget:
            for bx in stack:
                undecided.append({"box": [float(t) for t in bx],
                                  "area": float((bx[1]-bx[0])*(bx[3]-bx[2])),
                                  "reason": "budget"})
            break
        s0_, s1_, d0_, d1_ = stack.pop()
        T1 = iv.mpf([s0_, s1_])
        Dm = iv.mpf([d0_, d1_])
        T0 = T1 + Dm
        bb = {}
        ok = True
        for key, X in (("n0", T0), ("n1", T1), ("n0b", T0 - Tc.B),
                       ("n1b", T1 - Tc.B), ("nD", Dm)):
            ps = J.branch_pieces(X)
            if len(ps) != 1:
                ok = False
                break
            bb[key] = ps[0][1]
        br = bb if ok else None
        r = sep_res(Tc, T0, T1, br, Dm=Dm)
        if r is None:
            # massDet's ENCLOSURE straddles zero.  On D >= delta the true
            # massDet = pi^2 - phiJ(D)^2 is strictly positive (phiJ < pi for
            # D != 0 mod 2pi), so a straddle at a coarse box is pure interval
            # dependency: SUBDIVIDE, and record as undecided only at minwidth.
            if (s1_ - s0_) < minw and (d1_ - d0_) < minw:
                undecided.append({"box": [float(s0_), float(s1_), float(d0_),
                                          float(d1_)],
                                  "area": float((s1_-s0_)*(d1_-d0_)),
                                  "reason": "massDet_zero"})
                continue
            if s1_ - s0_ >= d1_ - d0_:
                m = (s0_ + s1_) / 2
                stack.append((s0_, m, d0_, d1_))
                stack.append((m, s1_, d0_, d1_))
            else:
                m = (d0_ + d1_) / 2
                stack.append((s0_, s1_, d0_, m))
                stack.append((s0_, s1_, m, d1_))
            continue
        if sgn(r[0]) != 0 or sgn(r[1]) != 0:
            continue
        k, N0, N1 = _krawczyk_sD(Tc, T1, Dm, br)
        if k == "empty":
            continue
        if k == "unique":
            S, DD = N0, N1
            for _ in range(60):
                kk, M0, M1 = _krawczyk_sD(Tc, S, DD, br)
                if kk != "unique" or (J.width(M0) >= J.width(S)
                                      and J.width(M1) >= J.width(DD)):
                    break
                S, DD = M0, M1
            R1 = S
            R0 = S + DD
            rr = sep_res(Tc, R0, R1, br, Dm=DD)
            lab, sg = schur_label(Tc, R0, R1, rr[2], rr[3])
            sols.append({
                "th0": [mpmath.nstr(J.lo(R0), 20), mpmath.nstr(J.hi(R0), 20)],
                "th1": [mpmath.nstr(J.lo(R1), 20), mpmath.nstr(J.hi(R1), 20)],
                "th0_mid": float(J.mid(R0)), "th1_mid": float(J.mid(R1)),
                "D_mid": float(J.mid(DD)),
                "enclosure_width": max(J.width(R1), J.width(DD)),
                "c0": [mpmath.nstr(J.lo(rr[2]), 12), mpmath.nstr(J.hi(rr[2]), 12)],
                "c1": [mpmath.nstr(J.lo(rr[3]), 12), mpmath.nstr(J.hi(rr[3]), 12)],
                "sgn_c0": sgn(rr[2]), "sgn_c1": sgn(rr[3]),
                "schur": lab, "schur_signs": list(sg),
            })
            continue
        if (s1_ - s0_) < minw and (d1_ - d0_) < minw:
            undecided.append({"box": [float(s0_), float(s1_), float(d0_),
                                      float(d1_)],
                              "area": float((s1_-s0_)*(d1_-d0_)),
                              "reason": "minwidth"})
            continue
        if s1_ - s0_ >= d1_ - d0_:
            m = (s0_ + s1_) / 2
            stack.append((s0_, m, d0_, d1_))
            stack.append((m, s1_, d0_, d1_))
        else:
            m = (d0_ + d1_) / 2
            stack.append((s0_, s1_, d0_, m))
            stack.append((s0_, s1_, m, d1_))
    collar = float(two_pi) * float(delta)
    return sols, undecided, steps, collar


def census_of(beta_f, y_f, sep=True, **kw):
    """The certified census string of the teacher (beta, massesAt(y))."""
    B = J.I(str(beta_f))
    Y = J.I(str(y_f))
    s0, s1 = masses_at(Y)
    Tc = Teacher(B, s0, s1)
    roots, undec = torque_roots(Tc)
    items = ["fit:global"]
    certified = not undec
    for r in roots:
        if r["label"] is None:
            certified = False
        elif r["label"]:
            items.append(r["label"])
    seps, sundec, steps, collar = ([], [], 0, 0.0)
    if sep:
        seps, sundec, steps, collar = separated_families(Tc, **kw)
        if sundec:
            certified = False
        # de-duplicate the s = 0 / s = 2pi seam (a family with th1 on the seam
        # is certified from both edge cells) and drop the EXACT FIT (both
        # students on the teacher lattice {0, beta}), which widgets.js lists
        # separately as fit:global and separatedScan drops.
        tp = 2 * math.pi
        bf = float(J.mid(Tc.B))

        def on_lattice(t):
            for cc in (0.0, bf % tp):
                d = (t - cc) % tp
                if min(d, tp - d) < 1e-6:
                    return True
            return False

        uniq = []
        for s in seps:
            key = (round(s["th1_mid"] % tp, 6), round(s["D_mid"], 6))
            if any(abs(key[0] - u) < 1e-5 and abs(key[1] - v) < 1e-5
                   for (u, v) in [(round(z["th1_mid"] % tp, 6),
                                   round(z["D_mid"], 6)) for z in uniq]):
                continue
            s["exact_fit"] = (on_lattice(s["th0_mid"])
                              and on_lattice(s["th1_mid"]))
            uniq.append(s)
        seps = uniq
        for s in seps:
            if s["exact_fit"]:
                continue
            if s["schur"] == "spurious":
                pos = s["sgn_c0"] >= 0 and s["sgn_c1"] >= 0
                # widgets.js: pos = c0 > -1e-9 and c1 > -1e-9
                items.append("separate:trap@positive" if pos
                             else "separate:trap@mixed")
            elif s["schur"].startswith("undecided_interval"):
                certified = False
    return {
        "beta": beta_f, "y": y_f,
        "n_torque_roots": len(roots),
        "torque_roots": roots,
        "torque_undecided": undec,
        "n_separated": len(seps),
        "separated": seps,
        "separated_undecided": sundec,
        "separated_steps": steps,
        "separated_collar_area": collar,
        "census": " | ".join(sorted(items)),
        "certified": certified,
    }


# --------------------------------------------------------------------------
# LOCATE-then-CERTIFY: a dense float locator (a GUESS, never a verdict) whose
# every candidate is then certified individually by Krawczyk, plus a separate
# exclusion pass that measures how much of the domain is certified solution
# free.  This separates the two obligations -- "every family listed is real"
# (certified per family) and "no family was missed" (certified by exclusion,
# with the uncovered area reported) -- instead of asking one branch and bound
# to do both.
# --------------------------------------------------------------------------

def _f_load(beta, s0, s1, t):
    return (s0 * J.f_phiJ(t) + s1 * J.f_phiJ(t - beta),
            -(s0 * J.f_HJ(t) + s1 * J.f_HJ(t - beta)))


def _f_res(beta, s0, s1, t0, t1):
    D = t0 - t1
    phiD = J.f_phiJ(D)
    hD = J.f_HJ(D)
    M = math.pi * math.pi - phiD * phiD
    if abs(M) < 1e-12:
        return None
    P0, A0 = _f_load(beta, s0, s1, t0)
    P1, A1 = _f_load(beta, s0, s1, t1)
    c0 = (math.pi * P0 - phiD * P1) / M
    c1 = (math.pi * P1 - phiD * P0) / M
    return (c1 * hD + A0, c0 * hD - A1, c0, c1)


def d_grid(dmin, nfine=70, ncoarse=170, dsplit=0.35):
    """D-values for the locator: GEOMETRIC from dmin to dsplit (so that
    near-coincident separated families, which is where this kernel puts its
    fine structure at beta near pi, are resolved at all), then uniform to pi."""
    out = []
    if dmin < dsplit:
        r = (dsplit / dmin) ** (1.0 / nfine)
        d = dmin
        for _ in range(nfine):
            out.append(d)
            d *= r
    a = max(dmin, dsplit)
    for k in range(ncoarse):
        out.append(a + (math.pi - a) * (k + 0.5) / ncoarse)
    return out


def locate_separated(beta, s0, s1, N=260, dmin=0.05, dvals=None):
    """Float grid + Newton polish over (th1, D), D in [dmin, pi].  A LOCATOR --
    a guess, never a verdict; every candidate is certified separately."""
    tp = 2 * math.pi
    found = []
    if dvals is None:
        dvals = d_grid(dmin)
    for i in range(N):
        for k in range(len(dvals)):
            t1 = tp * (i + 0.5) / N
            D = dvals[k]
            x = [t1 + D, t1]
            ok = False
            for _ in range(60):
                f = _f_res(beta, s0, s1, x[0], x[1])
                if f is None:
                    break
                if abs(f[0]) + abs(f[1]) < 1e-13:
                    ok = True
                    break
                h = 1e-7
                fa = _f_res(beta, s0, s1, x[0] + h, x[1])
                fb = _f_res(beta, s0, s1, x[0], x[1] + h)
                if fa is None or fb is None:
                    break
                j00 = (fa[0] - f[0]) / h
                j01 = (fb[0] - f[0]) / h
                j10 = (fa[1] - f[1]) / h
                j11 = (fb[1] - f[1]) / h
                det = j00 * j11 - j01 * j10
                if abs(det) < 1e-14:
                    break
                x = [x[0] + (-f[0] * j11 + f[1] * j01) / det,
                     x[1] + (-f[1] * j00 + f[0] * j10) / det]
            if not ok:
                continue
            a, b = x[0], x[1]
            d = a - b
            if d < 0:
                a, b, d = b, a, -d
            d = d % tp
            b = b % tp
            if d > math.pi:
                d = tp - d
                b = (b + tp - d) % tp
            if d < dmin:
                continue
            if any(abs(b - q[0]) < 1e-6 and abs(d - q[1]) < 1e-6 for q in found):
                continue
            found.append((b, d))
    found.sort()
    return found


def certify_candidate(Tc, s_f, D_f, r0=mpmath.mpf("1e-3")):
    """Krawczyk existence+uniqueness on a shrinking box around a located point."""
    r = r0
    for _ in range(40):
        S = iv.mpf([mpmath.mpf(repr(s_f)) - r, mpmath.mpf(repr(s_f)) + r])
        DD = iv.mpf([mpmath.mpf(repr(D_f)) - r, mpmath.mpf(repr(D_f)) + r])
        T0, T1 = S + DD, S
        bb = {}
        ok = True
        for key, X in (("n0", T0), ("n1", T1), ("n0b", T0 - Tc.B),
                       ("n1b", T1 - Tc.B), ("nD", DD)):
            ps = J.branch_pieces(X)
            if len(ps) != 1:
                ok = False
                break
            bb[key] = ps[0][1]
        br = bb if ok else None
        k, N0, N1 = _krawczyk_sD(Tc, S, DD, br)
        if k == "unique":
            for _ in range(80):
                kk, M0, M1 = _krawczyk_sD(Tc, N0, N1, br)
                if kk != "unique":
                    break
                if J.width(M0) >= J.width(N0) and J.width(M1) >= J.width(N1):
                    break
                N0, N1 = M0, M1
            R1, DDf = N0, N1
            R0 = R1 + DDf
            rr = sep_res(Tc, R0, R1, br, Dm=DDf)
            lab, sg = schur_label(Tc, R0, R1, rr[2], rr[3])
            return {
                "th0": [mpmath.nstr(J.lo(R0), 22), mpmath.nstr(J.hi(R0), 22)],
                "th1": [mpmath.nstr(J.lo(R1), 22), mpmath.nstr(J.hi(R1), 22)],
                "th0_mid": float(J.mid(R0)), "th1_mid": float(J.mid(R1)),
                "D_mid": float(J.mid(DDf)),
                "enclosure_width": max(J.width(R1), J.width(DDf)),
                "c0": [mpmath.nstr(J.lo(rr[2]), 14), mpmath.nstr(J.hi(rr[2]), 14)],
                "c1": [mpmath.nstr(J.lo(rr[3]), 14), mpmath.nstr(J.hi(rr[3]), 14)],
                "sgn_c0": sgn(rr[2]), "sgn_c1": sgn(rr[3]),
                "schur": lab, "schur_signs": list(sg),
                "krawczyk_box_radius": float(r),
            }
        r = r / 2
    return None


def separated_certified(beta_f, y_f, N=260, dmin=1e-4, dvals=None):
    """Located-then-certified separated families of one teacher."""
    B = J.I(str(beta_f))
    Y = J.I(str(y_f))
    s0, s1 = masses_at(Y)
    Tc = Teacher(B, s0, s1)
    f0, f1 = float(J.mid(s0)), float(J.mid(s1))
    cands = locate_separated(beta_f, f0, f1, N=N, dmin=dmin, dvals=dvals)
    out = []
    failed = []
    for (s_, d_) in cands:
        rec = certify_candidate(Tc, s_, d_)
        if rec is None:
            failed.append((s_, d_))
        else:
            out.append(rec)
    # de-duplicate on the UNORDERED pair modulo 2pi, and drop the EXACT FIT
    # (both students on the teacher lattice {0, beta} mod 2pi), which
    # widgets.js lists separately as fit:global and separatedScan drops.
    tp = 2 * math.pi

    def key(rec):
        a, b = rec["th0_mid"] % tp, rec["th1_mid"] % tp
        return tuple(sorted((round(a, 7), round(b, 7))))

    def on_lattice(t):
        for c in (0.0, beta_f % tp):
            d = (t - c) % tp
            if min(d, tp - d) < 1e-5:
                return True
        return False

    uniq, dropped_fit = [], []
    for rec in out:
        k = key(rec)
        if any(abs(k[0] - u[0]) < 1e-5 and abs(k[1] - u[1]) < 1e-5
               for u in [key(z) for z in uniq] + [key(z) for z in dropped_fit]):
            continue
        if on_lattice(rec["th0_mid"]) and on_lattice(rec["th1_mid"]):
            rec["exact_fit"] = True
            dropped_fit.append(rec)
            continue
        uniq.append(rec)
    return Tc, uniq, failed, len(cands), dropped_fit


def certified_census(beta_f, y_f, N=200, dmin=1e-4, tol_root=None):
    """The full certified census string of one teacher."""
    Tc, seps, failed, ncand, fits = separated_certified(beta_f, y_f, N=N,
                                                        dmin=dmin)
    seps = [z for z in seps if z["D_mid"] > 1e-4]
    roots, undec = torque_roots(Tc)
    items = ["fit:global"]
    certified = not undec and not failed
    for r in roots:
        if r["label"] is None:
            certified = False
        elif r["label"]:
            items.append(r["label"])
    for s_ in seps:
        if s_["schur"] == "spurious":
            pos = s_["sgn_c0"] >= 0 and s_["sgn_c1"] >= 0
            items.append("separate:trap@positive" if pos
                         else "separate:trap@mixed")
        elif s_["schur"].startswith("undecided"):
            certified = False
    return {
        "beta": beta_f, "y": y_f,
        "n_torque_roots": len(roots),
        "torque_roots": roots,
        "torque_undecided": undec,
        "n_separated": len(seps),
        "separated": seps,
        "n_exact_fit_dropped": len(fits),
        "locator_candidates": ncand,
        "locator_failed_certification": failed,
        "census": " | ".join(sorted(items)),
        "certified_modulo_locator_completeness": certified,
    }
