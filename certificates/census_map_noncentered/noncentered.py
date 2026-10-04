"""
NONCENTERED (plain-ReLU) two-student kernel: float mirror of
website/widgets.js kernelOf("noncentered") + rigorous interval arithmetic
(mpmath.iv) implementations.

Working precision is set by set_prec(); everything certified goes through
mpmath.iv (directed rounding).  float64 is used ONLY in the float mirror,
which exists to be compared bit-for-bit against widgets.js, and in pre-scans
that GUESS where roots are (never in a reported enclosure).

Kernel (widgets.js lines 231/232/233/294, kernelOf line 238):
    phiJ(t) = (pi-x) cos x + sin x,   x = mod(t,2pi) folded by x>pi -> 2pi-x
    HJ(t)   = (pi-x) sin x  (x<=pi) ,  (x-pi) sin x  (x>pi),  x = mod(t,2pi)
    dHJ(t)  = HJ'(t) = phiJ(t) - 2|sin t| = slopeAtom
    kap = pi,  per = 2pi
"""
import math
import mpmath
from mpmath import iv, mp

PREC_BITS = 200          # working precision of the certified path


def set_prec(bits=PREC_BITS):
    iv.prec = bits
    mp.prec = bits + 20


set_prec()

PI = math.pi

# --------------------------------------------------------------------------
# float mirror of widgets.js -- NOT part of the certified path.
# Arithmetic order copied verbatim from the source lines.
# --------------------------------------------------------------------------


def f_mod(x, m):
    r = math.fmod(x, m)
    return r + m if r < 0 else r


def f_phiJ(t):
    x = f_mod(t, 2 * PI)
    if x > PI:
        x = 2 * PI - x
    return (PI - x) * math.cos(x) + math.sin(x)


def f_HJ(t):
    x = f_mod(t, 2 * PI)
    return (PI - x) * math.sin(x) if x <= PI else (x - PI) * math.sin(x)


def f_dHJ(t):
    x = f_mod(t, 2 * PI)
    return ((PI - x) * math.cos(x) - math.sin(x) if x <= PI
            else math.sin(x) + (x - PI) * math.cos(x))


def f_slopeAtom(t):
    return f_phiJ(t) - 2 * abs(math.sin(t))


def f_Wtau(beta, t):
    return f_HJ(t) * f_slopeAtom(t - beta) - f_HJ(t - beta) * f_slopeAtom(t)


def f_Wpot(beta, t):
    return f_phiJ(t) * f_HJ(t - beta) - f_phiJ(t - beta) * f_HJ(t)


def f_Wwgt(beta, t):
    return abs(math.sin(t)) * f_HJ(t - beta) - abs(math.sin(t - beta)) * f_HJ(t)


KAP = PI


def f_balances(beta, s0, s1, t0, t1):
    """widgets.js separatedScan residual, h(D) guard removed."""
    D = t0 - t1
    H = f_HJ(D)
    def Trq(t):
        return -(s0 * f_HJ(t) + s1 * f_HJ(t - beta))
    def Pot(t):
        return s0 * f_phiJ(t) + s1 * f_phiJ(t - beta)
    c0 = Trq(t1) / H
    c1 = -Trq(t0) / H
    return (KAP * c0 + f_phiJ(D) * c1 - Pot(t0),
            f_phiJ(D) * c0 + KAP * c1 - Pot(t1))


def f_Dsep(beta, t0, t1):
    """H(D)^2 * separatedDet, in the division form (matches ref_widgetsJ.js)."""
    A = f_balances(beta, 1.0, 0.0, t0, t1)
    B = f_balances(beta, 0.0, 1.0, t0, t1)
    HD = f_HJ(t0 - t1)
    return (A[0] * B[1] - B[0] * A[1]) * HD * HD


def f_ratioCoord(s0, s1):
    a = math.atan2(s0, s1)
    a = f_mod(a, PI)
    if a <= 0:
        a += PI
    return 2 * a / PI - 1


# --------------------------------------------------------------------------
# interval helpers
# --------------------------------------------------------------------------

def I(a, b=None):
    if b is None:
        b = a
    return iv.mpf([a, b])


def PI_IV():
    return iv.pi


def endpoints(z):
    """Exact mp.mpf endpoints of an interval."""
    a = mpmath.mpf()
    a._mpf_ = z._mpi_[0]
    b = mpmath.mpf()
    b._mpf_ = z._mpi_[1]
    return a, b


def hull(*xs):
    if len(xs) == 1:
        return xs[0]
    los, his = [], []
    for x in xs:
        a, b = endpoints(x)
        los.append(a)
        his.append(b)
    return iv.mpf([min(los), max(his)])


def contains_zero(x):
    a, b = endpoints(x)
    return a <= 0 <= b


def width(x):
    a, b = endpoints(x)
    return float(b - a)


def mid(x):
    a, b = endpoints(x)
    return (a + b) / 2


def lo(x):
    return endpoints(x)[0]


def hi(x):
    return endpoints(x)[1]


def rad(x):
    a, b = endpoints(x)
    return float(max(abs(a), abs(b)))


# --------------------------------------------------------------------------
# HALF-BRANCH atoms.
#
# Half-branch n:  t in [n*pi, (n+1)*pi].  Put
#     p = pi + 2pi*floor(n/2) - t          (so x = mod(t,2pi) = pi - p)
#     e = +1 if n even else -1             (e = sign(p); p in [0,pi] / [-pi,0])
#     P = e*p = |p| ,  s = sin p , c = cos p
# Then the CLOSED BRANCH FORMS of widgets.js are exactly
#     h   = P*s
#     phi = -P*c + e*s
#     sA  = -P*c - e*s   ( = phi - 2|sin t| = h' )
#     |sin t| = e*s
# and, since dp/dt = -1,
#     h'   = sA
#     phi' = -h
#     sA'  = 2*e*c - P*s
#     (|sin t|)' = -e*c
# These are ANALYTIC on the closed half-branch; the kinks of the kernel are
# exactly the half-branch endpoints p in {0, +-pi}.
# --------------------------------------------------------------------------

class Atoms(object):
    __slots__ = ("phi", "h", "sA", "asin", "dphi", "dh", "dsA", "dasin", "p")


def p_of(T, n):
    """p-coordinate of t on half-branch n (T an interval)."""
    return PI_IV() * (1 + 2 * (n // 2)) - T


def atoms_on_branch(T, n):
    """T: interval of t known to lie in [n*pi, (n+1)*pi] (small overshoot ok)."""
    p = p_of(T, n)
    e = 1 if n % 2 == 0 else -1
    P = e * p
    s = iv.sin(p)
    c = iv.cos(p)
    A = Atoms()
    A.p = p
    A.h = P * s
    A.phi = -(P * c) + e * s
    A.sA = -(P * c) - e * s
    A.asin = e * s
    A.dh = A.sA
    A.dphi = -A.h
    A.dsA = 2 * e * c - P * s
    A.dasin = -(e * c)
    return A


def branch_pieces(T):
    """Rigorous cover of interval T by (piece, n) with piece inside half-branch n."""
    Ta, Tb = endpoints(T)
    nlo = int(mpmath.floor(Ta / mpmath.mpf(3.14159265358979))) - 2
    nhi = int(mpmath.floor(Tb / mpmath.mpf(3.14159265358980))) + 2
    out = []
    for n in range(nlo, nhi + 1):
        a = endpoints(iv.mpf(n) * PI_IV())[0]          # <= n*pi
        b = endpoints(iv.mpf(n + 1) * PI_IV())[1]      # >= (n+1)*pi
        plo = max(Ta, a)
        phi_ = min(Tb, b)
        if plo <= phi_:
            out.append((iv.mpf([plo, phi_]), n))
    assert out, "empty branch cover"
    assert endpoints(out[0][0])[0] == Ta and endpoints(out[-1][0])[1] == Tb
    for i in range(len(out) - 1):
        assert endpoints(out[i + 1][0])[0] <= endpoints(out[i][0])[1]
    return out


def atoms(T):
    """Enclosure of all atoms over T, valid across kinks (hull over branches).

    Returns (Atoms, smooth_flag).  smooth_flag is True iff T met exactly one
    half-branch, i.e. no kink was crossed."""
    ps = branch_pieces(T)
    As = [atoms_on_branch(p, n) for p, n in ps]
    R = Atoms()
    for f in Atoms.__slots__:
        setattr(R, f, hull(*[getattr(A, f) for A in As]))
    return R, (len(ps) == 1)


# --------------------------------------------------------------------------
# The three mass-free determinants.
# --------------------------------------------------------------------------

def _dets(At, Au):
    Wt = At.h * Au.sA - Au.h * At.sA
    Wp = At.phi * Au.h - Au.phi * At.h
    Ww = At.asin * Au.h - Au.asin * At.h
    dWt = At.h * Au.dsA - Au.h * At.dsA
    dWp = At.phi * Au.sA - Au.phi * At.sA
    dWw = At.dasin * Au.h + At.asin * Au.dh - Au.dasin * At.h - Au.asin * At.dh
    return Wt, Wp, Ww, dWt, dWp, dWw


def dets_branch(B, T, nt, nu):
    """Branch indices supplied by the caller (legitimate when T is inside the
    CLOSED half-branch nt and T-B inside nu)."""
    At = atoms_on_branch(T, nt)
    Au = atoms_on_branch(T - B, nu)
    Wt, Wp, Ww, dWt, dWp, dWw = _dets(At, Au)
    return {"Wtau": (Wt, dWt), "Wpot": (Wp, dWp), "Wwgt": (Ww, dWw)}


def dets_iv(B, T):
    At, st = atoms(T)
    Au, su = atoms(T - B)
    Wt, Wp, Ww, dWt, dWp, dWw = _dets(At, Au)
    sm = st and su
    return {"Wtau": (Wt, dWt, sm), "Wpot": (Wp, dWp, sm), "Wwgt": (Ww, dWw, sm)}


# --------------------------------------------------------------------------
# REDUCED (kink-free) form.
#
# On a piece where t sits in one half-branch and t-beta in one half-branch,
# with p = p(t), q = p(t-beta), d = q - p (an element of beta + 2pi Z) and
# E = e_t * e_u = sign(p)*sign(q),
#
#     Wtau = E * F(p,q),  Wpot = E * G(p,q),  Wwgt = E * N(p,q)
#     Qd(p,q) = p q sin(q-p)
#     N (p,q) = (q-p) sin p sin q
#     F = Qd + N ,  G = -Qd + N          [hence F + G = 2N identically]
#
# In the midpoint variable m = (p+q)/2 = p + d/2  (so p = m-d/2, q = m+d/2):
#     Qd = (m^2 - d^2/4) sin d
#     N  = (d/2)(cos d - cos 2m)
# both ENTIRE in m: the reduction absorbed every kink into the piece boundary.
# --------------------------------------------------------------------------

def redu(d, m):
    """d, m intervals -> (F, G, N, F_m, G_m, N_m)."""
    sd = iv.sin(d)
    cd = iv.cos(d)
    Qd = (m * m - d * d / 4) * sd
    Nn = (d / 2) * (cd - iv.cos(2 * m))
    Qm = 2 * m * sd
    Nm = d * iv.sin(2 * m)
    return (Qd + Nn, -Qd + Nn, Nn, Qm + Nm, -Qm + Nm, Nm)


def redu_mm(d, m):
    """second m-derivatives (F_mm, G_mm, N_mm)."""
    sd = iv.sin(d)
    Qmm = 2 * sd
    Nmm = 2 * d * iv.cos(2 * m)
    return (Qmm + Nmm, -Qmm + Nmm, Nmm)


def pieces_of(beta_iv):
    """The four smooth pieces of the t-period [0,2pi) for beta in (0,2pi),
    beta != pi, described in the (d, m) variables.

    Returns a list of dicts with
       'name', 'd' (interval), 'm_lo', 'm_hi' (intervals), 'E' (+-1),
       'nt', 'nu' (half-branch indices of t and t-beta),
       't_lo', 't_hi' (interval endpoints of the t-piece)
    The pieces tile [0, 2pi] in t, in increasing t order.

    p = pi - t on [0,2pi];  q = p + d;  m = p + d/2.
    Branch d = beta      on p in [-pi, pi-beta]   (t in [beta, 2pi])
    Branch d = beta-2pi  on p in [pi-beta, pi]    (t in [0, beta])
    each split further at p = 0 and q = 0.
    """
    PIv = PI_IV()
    b = beta_iv
    out = []
    # kink t-values: 0, beta, pi, beta+pi (mod 2pi), 2pi
    kt = [I(0), b, PIv, b + PIv if hi(b) < hi(PIv) else b - PIv, 2 * PIv]
    kt = sorted(kt, key=lambda z: float(mid(z)))
    for i in range(len(kt) - 1):
        tlo, thi = kt[i], kt[i + 1]
        tm = (tlo + thi) / 2
        tmf = float(mid(tm))
        # half-branch of t and of t-beta at the midpoint
        nt = int(math.floor(tmf / PI))
        nu = int(math.floor((tmf - float(mid(b))) / PI))
        # p-range: p = pi + 2pi*floor(nt/2) - t
        p_hi = p_of(tlo, nt)
        p_lo = p_of(thi, nt)
        q_hi = p_of(tlo - b, nu)
        q_lo = p_of(thi - b, nu)
        d = (q_lo - p_lo)          # constant on the piece
        m_lo = (p_lo + q_lo) / 2
        m_hi = (p_hi + q_hi) / 2
        e_t = 1 if nt % 2 == 0 else -1
        e_u = 1 if nu % 2 == 0 else -1
        out.append({"t_lo": tlo, "t_hi": thi, "nt": nt, "nu": nu,
                    "d": d, "m_lo": m_lo, "m_hi": m_hi, "E": e_t * e_u,
                    "p_lo": p_lo, "p_hi": p_hi, "q_lo": q_lo, "q_hi": q_hi})
    return out


def t_from_m(pc, m):
    """t from the midpoint variable m on piece pc."""
    # m = p + d/2 and p = pi + 2pi*floor(nt/2) - t
    return PI_IV() * (1 + 2 * (pc["nt"] // 2)) - (m - pc["d"] / 2)


# --------------------------------------------------------------------------
# mass direction / y coordinate
#   (s0, s1) proportional to (-h(t-beta), h(t)),  y = ratioCoord(s0,s1)
# --------------------------------------------------------------------------

def y_from_s(s0, s1):
    """Certified y enclosure from an (s0,s1) enclosure (widgets.js ratioCoord).
    Returns None if the enclosure straddles the fold seam."""
    if contains_zero(s0) and contains_zero(s1):
        return None
    a = iv.atan2(s0, s1)
    if lo(a) > 0:
        psi = a
    elif hi(a) < 0:
        psi = a + PI_IV()
    else:
        return None
    return 2 * psi / PI_IV() - 1


# --------------------------------------------------------------------------
# SEPARATED stratum.
#
# D = th0 - th1, hD = h(D), phiD = phi(D), kap = pi.  With
#   c0 = A(th1)/hD, c1 = -A(th0)/hD,
#   F0 = kap c0 + phiD c1 - P(th0),   F1 = phiD c0 + kap c1 - P(th1),
# the columns against the unit teachers e0=(1,0), e1=(0,1) are (times hD)
#   a0 = -kap h(th1) + phiD h(th0) - hD phi(th0)      [= hD * F0(e0)]
#   b0 = -phiD h(th1) + kap h(th0) - hD phi(th1)      [= hD * F1(e0)]
#   a1 = -kap h(th1-B) + phiD h(th0-B) - hD phi(th0-B)
#   b1 = -phiD h(th1-B) + kap h(th0-B) - hD phi(th1-B)
# and
#   Dsep := a0*b1 - a1*b0 = hD^2 * separatedDet .
# Dsep involves only h and phi, which are C^1 everywhere, so Dsep is C^1.
# Kernel vector (mass direction): (s0,s1) ~ (-a1, a0).
# --------------------------------------------------------------------------

class Sep(object):
    __slots__ = ("D", "d0", "d1", "a0", "b0", "a1", "b1")


def sep_atoms(B, T0, T1, branches=None):
    """branches: optional dict with keys 'n0','n1','n0b','n1b','nD' giving the
    half-branch indices of th0, th1, th0-B, th1-B, th0-th1.  If given, the
    closed branch forms are used (legitimate on the closed half-branch);
    otherwise the multi-branch hull is used."""
    Dm = T0 - T1
    if branches is None:
        A0, s1_ = atoms(T0)
        A1, s2_ = atoms(T1)
        A0b, s3_ = atoms(T0 - B)
        A1b, s4_ = atoms(T1 - B)
        AD, s5_ = atoms(Dm)
        smooth = s1_ and s2_ and s3_ and s4_ and s5_
    else:
        A0 = atoms_on_branch(T0, branches["n0"])
        A1 = atoms_on_branch(T1, branches["n1"])
        A0b = atoms_on_branch(T0 - B, branches["n0b"])
        A1b = atoms_on_branch(T1 - B, branches["n1b"])
        AD = atoms_on_branch(Dm, branches["nD"])
        smooth = True
    return A0, A1, A0b, A1b, AD, smooth


def sep_det(B, T0, T1, branches=None):
    """Returns (Dsep, a0, a1, smooth)."""
    A0, A1, A0b, A1b, AD, smooth = sep_atoms(B, T0, T1, branches)
    kap = PI_IV()
    hD, phiD = AD.h, AD.phi
    a0 = -(kap * A1.h) + phiD * A0.h - hD * A0.phi
    b0 = -(phiD * A1.h) + kap * A0.h - hD * A1.phi
    a1 = -(kap * A1b.h) + phiD * A0b.h - hD * A0b.phi
    b1 = -(phiD * A1b.h) + kap * A0b.h - hD * A1b.phi
    return a0 * b1 - a1 * b0, a0, a1, smooth


def sep_det_grad(B, T0, T1, branches=None):
    """Returns (Dsep, dDsep/dth0, dDsep/dth1, a0, a1, smooth)."""
    A0, A1, A0b, A1b, AD, smooth = sep_atoms(B, T0, T1, branches)
    kap = PI_IV()
    hD, phiD = AD.h, AD.phi
    sAD = AD.sA                          # h'(D)
    a0 = -(kap * A1.h) + phiD * A0.h - hD * A0.phi
    b0 = -(phiD * A1.h) + kap * A0.h - hD * A1.phi
    a1 = -(kap * A1b.h) + phiD * A0b.h - hD * A0b.phi
    b1 = -(phiD * A1b.h) + kap * A0b.h - hD * A1b.phi
    # First partials, with the two cancelling terms removed ALGEBRAICALLY
    # (phi'(D) = -h(D) makes dphiD*h0 and -hD*dphi0 cancel exactly).  Keeping
    # them would leave a dependency blow-up of one to two orders of magnitude
    # in the interval enclosure of the gradient.
    #   a0_0 =  phiD*sA0  - sAD*phi0
    #   b0_1 = -phiD*sA1  + sAD*phi1
    # and likewise for the beta-shifted column.
    a0_0 = phiD * A0.sA - sAD * A0.phi
    b0_1 = -(phiD * A1.sA) + sAD * A1.phi
    a1_0 = phiD * A0b.sA - sAD * A0b.phi
    b1_1 = -(phiD * A1b.sA) + sAD * A1b.phi
    a0_1 = -(kap * A1.sA) + hD * A0.h + sAD * A0.phi
    b0_0 = hD * A1.h + kap * A0.sA - sAD * A1.phi
    a1_1 = -(kap * A1b.sA) + hD * A0b.h + sAD * A0b.phi
    b1_0 = hD * A1b.h + kap * A0b.sA - sAD * A1b.phi
    # d/dbeta : only the "-B" atoms move, with d/dbeta = -d/dt
    a1_b = kap * A1b.sA - phiD * A0b.sA - hD * A0b.h
    b1_b = phiD * A1b.sA - kap * A0b.sA - hD * A1b.h
    Dv = a0 * b1 - a1 * b0
    D0 = a0_0 * b1 + a0 * b1_0 - a1_0 * b0 - a1 * b0_0
    D1 = a0_1 * b1 + a0 * b1_1 - a1_1 * b0 - a1 * b0_1
    Db = a0 * b1_b - a1_b * b0
    return Dv, D0, D1, Db, a0, a1, AD, smooth


def sep_columns(B, T0, T1, branches=None):
    """All four column entries and all their first partials, plus the D-atoms."""
    A0, A1, A0b, A1b, AD, smooth = sep_atoms(B, T0, T1, branches)
    kap = PI_IV()
    hD, phiD, sAD = AD.h, AD.phi, AD.sA
    C = {}
    C["a0"] = -(kap * A1.h) + phiD * A0.h - hD * A0.phi
    C["b0"] = -(phiD * A1.h) + kap * A0.h - hD * A1.phi
    C["a1"] = -(kap * A1b.h) + phiD * A0b.h - hD * A0b.phi
    C["b1"] = -(phiD * A1b.h) + kap * A0b.h - hD * A1b.phi
    C["a0_0"] = phiD * A0.sA - sAD * A0.phi
    C["b0_1"] = -(phiD * A1.sA) + sAD * A1.phi
    C["a1_0"] = phiD * A0b.sA - sAD * A0b.phi
    C["b1_1"] = -(phiD * A1b.sA) + sAD * A1b.phi
    C["a0_1"] = -(kap * A1.sA) + hD * A0.h + sAD * A0.phi
    C["b0_0"] = hD * A1.h + kap * A0.sA - sAD * A1.phi
    C["a1_1"] = -(kap * A1b.sA) + hD * A0b.h + sAD * A0b.phi
    C["b1_0"] = hD * A1b.h + kap * A0b.sA - sAD * A1b.phi
    C["a0_b"] = iv.mpf(0)
    C["b0_b"] = iv.mpf(0)
    C["a1_b"] = kap * A1b.sA - phiD * A0b.sA - hD * A0b.h
    C["b1_b"] = phiD * A1b.sA - kap * A0b.sA - hD * A1b.h
    return C, AD, smooth


def _inter(x, y):
    a1, b1 = endpoints(x)
    a2, b2 = endpoints(y)
    a, b = max(a1, a2), min(b1, b2)
    if a > b:
        return x
    return iv.mpf([a, b])


# --------------------------------------------------------------------------
# WELL-CONDITIONED mass-free determinant (the radial-row elimination).
#
# Solve the student weights from the RADIAL rows (massDet = kap^2 - phi(D)^2,
# which vanishes ONLY at D = 0 mod 2pi), substitute into the two TORQUE rows,
# and clear the division by massDet.  Against the unit teachers e0 = (1,0),
# e1 = (0,1) the two cleared torque residuals are (M = massDet, kap = pi)
#
#   a0 = (kap phi(th1) - phi(D) phi(th0)) h(D) - M h(th0)      [ = M G0(e0) ]
#   b0 = (kap phi(th0) - phi(D) phi(th1)) h(D) + M h(th1)      [ = M G1(e0) ]
#   a1, b1 : the same with every atom argument shifted by -beta [ = M G(e1) ]
#
#   DsepWC := a0 b1 - a1 b0.
#
# A separated critical configuration of a teacher (s0,s1) != (0,0) at
# (th0, th1) with massDet != 0 forces DsepWC = 0 (the same linearity-in-the-
# masses argument as separatedDet_eq_zero_of_isSeparatedBalanced, run on the
# radial-solved residual instead of the torque-solved one), and conversely a
# zero of DsepWC at which the two linear forms are not both identically zero
# is realised by the kernel-vector teacher (s0,s1) ~ (-a1, a0).
# DsepWC is NOT proportional to Dsep: the two eliminations acquire different
# extraneous loci (Dsep degenerates on h(D) = 0, i.e. D = 0 mod pi; DsepWC on
# massDet = 0, i.e. D = 0 mod 2pi only).  On the separated stratum both vanish.
# --------------------------------------------------------------------------

def wc_columns(B, T0, T1, branches=None):
    """The four WC column entries and their first partials, plus the D-atoms."""
    A0, A1, A0b, A1b, AD, smooth = sep_atoms(B, T0, T1, branches)
    kap = PI_IV()
    hD, phiD, sAD = AD.h, AD.phi, AD.sA
    M = kap * kap - phiD * phiD
    # dM/dth0 = -2 phiD dphiD/dth0 = -2 phiD (-hD) = 2 phiD hD ; dM/dth1 = -that
    Mp = 2 * phiD * hD
    C = {}
    C["a0"] = (kap * A1.phi - phiD * A0.phi) * hD - M * A0.h
    C["b0"] = (kap * A0.phi - phiD * A1.phi) * hD + M * A1.h
    C["a1"] = (kap * A1b.phi - phiD * A0b.phi) * hD - M * A0b.h
    C["b1"] = (kap * A0b.phi - phiD * A1b.phi) * hD + M * A1b.h
    # partials.  dphi = -h, dh = sA; D-atom derivatives w.r.t. th0 carry +1,
    # w.r.t. th1 carry -1; the "-beta" atoms move with beta by -d/dt.
    C["a0_0"] = ((hD * A0.phi + phiD * A0.h) * hD
                 + (kap * A1.phi - phiD * A0.phi) * sAD
                 - Mp * A0.h - M * A0.sA)
    C["a0_1"] = (-(kap * A1.h) - hD * A0.phi) * hD \
        - (kap * A1.phi - phiD * A0.phi) * sAD + Mp * A0.h
    C["b0_0"] = (-(kap * A0.h) + hD * A1.phi) * hD \
        + (kap * A0.phi - phiD * A1.phi) * sAD + Mp * A1.h
    C["b0_1"] = ((phiD * A1.h - hD * A1.phi) * hD
                 - (kap * A0.phi - phiD * A1.phi) * sAD
                 - Mp * A1.h + M * A1.sA)
    C["a1_0"] = ((hD * A0b.phi + phiD * A0b.h) * hD
                 + (kap * A1b.phi - phiD * A0b.phi) * sAD
                 - Mp * A0b.h - M * A0b.sA)
    C["a1_1"] = (-(kap * A1b.h) - hD * A0b.phi) * hD \
        - (kap * A1b.phi - phiD * A0b.phi) * sAD + Mp * A0b.h
    C["b1_0"] = (-(kap * A0b.h) + hD * A1b.phi) * hD \
        + (kap * A0b.phi - phiD * A1b.phi) * sAD + Mp * A1b.h
    C["b1_1"] = ((phiD * A1b.h - hD * A1b.phi) * hD
                 - (kap * A0b.phi - phiD * A1b.phi) * sAD
                 - Mp * A1b.h + M * A1b.sA)
    # d/dbeta: only the "-beta" atoms move, with d/dbeta = -d/dt
    C["a0_b"] = iv.mpf(0)
    C["b0_b"] = iv.mpf(0)
    C["a1_b"] = (kap * A1b.h - phiD * A0b.h) * hD + M * A0b.sA
    C["b1_b"] = (kap * A0b.h - phiD * A1b.h) * hD - M * A1b.sA
    return C, AD, M, smooth


def wc_det_grad_mv(B, T0, T1, branches=None):
    """DsepWC, its gradient and the kernel-vector pair (a0, a1), with the four
    column entries MEAN-VALUE refined against their own first partials.  Sound
    because the columns are C^1 in (beta, th0, th1) (h and phi are C^1).

    Returns (Dv, D0, D1, Db, a0, a1, AD, M, smooth, Dc)."""
    C, AD, M, smooth = wc_columns(B, T0, T1, branches)
    bm = iv.mpf([mid(B), mid(B)])
    t0m = iv.mpf([mid(T0), mid(T0)])
    t1m = iv.mpf([mid(T1), mid(T1)])
    Cc, ADc, Mc, _ = wc_columns(bm, t0m, t1m, branches)
    rb, r0, r1 = B - bm, T0 - t0m, T1 - t1m
    V = {}
    for nm in ("a0", "b0", "a1", "b1"):
        V[nm] = _inter(C[nm], Cc[nm] + C[nm + "_0"] * r0 + C[nm + "_1"] * r1
                       + C[nm + "_b"] * rb)
    a0, b0, a1, b1 = V["a0"], V["b0"], V["a1"], V["b1"]
    Dv = a0 * b1 - a1 * b0
    D0 = C["a0_0"] * b1 + a0 * C["b1_0"] - C["a1_0"] * b0 - a1 * C["b0_0"]
    D1 = C["a0_1"] * b1 + a0 * C["b1_1"] - C["a1_1"] * b0 - a1 * C["b0_1"]
    Db = a0 * C["b1_b"] - C["a1_b"] * b0 - a1 * C["b0_b"]
    Dc = Cc["a0"] * Cc["b1"] - Cc["a1"] * Cc["b0"]
    Dv = _inter(Dv, Dc + D0 * r0 + D1 * r1 + Db * rb)
    return Dv, D0, D1, Db, a0, a1, AD, M, smooth, Dc


def f_DsepWC(beta, t0, t1):
    """float mirror of DsepWC (locator/diagnostic only, never a verdict)."""
    kap = PI
    D = t0 - t1
    phiD, hD = f_phiJ(D), f_HJ(D)
    M = kap * kap - phiD * phiD
    a0 = (kap * f_phiJ(t1) - phiD * f_phiJ(t0)) * hD - M * f_HJ(t0)
    b0 = (kap * f_phiJ(t0) - phiD * f_phiJ(t1)) * hD + M * f_HJ(t1)
    a1 = (kap * f_phiJ(t1 - beta) - phiD * f_phiJ(t0 - beta)) * hD \
        - M * f_HJ(t0 - beta)
    b1 = (kap * f_phiJ(t0 - beta) - phiD * f_phiJ(t1 - beta)) * hD \
        + M * f_HJ(t1 - beta)
    return a0 * b1 - a1 * b0


# --------------------------------------------------------------------------
# REDUCED well-conditioned determinant.  sympy-exact identity
# (validate_wc.py):
#
#   DsepWC = massDet * Dred,
#   Dred   = hD^2 Wphi + hD (kap (Wp0 + Wp1) - phiD X) - massDet Wh
#
# with the five two-atom WRONSKIAN GROUPS (0 = th0, 1 = th1, b = "-beta")
#   Wphi = phi1 phi0b - phi0 phi1b
#   Wp0  = phi0 h0b - phi0b h0          [ = potentialAtomWronskian(beta, th0) ]
#   Wp1  = phi1 h1b - phi1b h1          [ = potentialAtomWronskian(beta, th1) ]
#   X    = phi0 h1b + phi1 h0b - phi1b h0 - phi0b h1
#   Wh   = h0 h1b - h0b h1
#
# Each group's first partials enjoy EXACT cancellations (phi' = -h, h' = sA
# kill the h.h cross terms of the potential Wronskians), so a mean-value
# refinement at group level loses far less to dependency than the same
# refinement applied to the raw 4-column product.  massDet > 0 off the
# 2pi-coincidence lattice, so {Dred = 0} = {DsepWC = 0} there.
# --------------------------------------------------------------------------

def _wc_groups(B, T0, T1, branches=None):
    """The five groups, D-atoms, and all group partials (th0, th1, beta)."""
    A0, A1, A0b, A1b, AD, smooth = sep_atoms(B, T0, T1, branches)
    G = {}
    G["Wphi"] = A1.phi * A0b.phi - A0.phi * A1b.phi
    G["Wp0"] = A0.phi * A0b.h - A0b.phi * A0.h
    G["Wp1"] = A1.phi * A1b.h - A1b.phi * A1.h
    G["X"] = A0.phi * A1b.h + A1.phi * A0b.h - A1b.phi * A0.h - A0b.phi * A1.h
    G["Wh"] = A0.h * A1b.h - A0b.h * A1.h
    # partials, with the exact cancellations applied
    G["Wphi_0"] = A0.h * A1b.phi - A0b.h * A1.phi
    G["Wphi_1"] = A0.phi * A1b.h - A0b.phi * A1.h
    G["Wphi_b"] = A1.phi * A0b.h - A0.phi * A1b.h
    G["Wp0_0"] = A0.phi * A0b.sA - A0b.phi * A0.sA
    G["Wp0_1"] = iv.mpf(0)
    G["Wp0_b"] = -(A0.phi * A0b.sA) - A0b.h * A0.h
    G["Wp1_0"] = iv.mpf(0)
    G["Wp1_1"] = A1.phi * A1b.sA - A1b.phi * A1.sA
    G["Wp1_b"] = -(A1.phi * A1b.sA) - A1b.h * A1.h
    G["X_0"] = (A0b.h * A1.h - A0.h * A1b.h) \
        + (A1.phi * A0b.sA - A1b.phi * A0.sA)
    G["X_1"] = (A0.h * A1b.h - A0b.h * A1.h) \
        + (A0.phi * A1b.sA - A0b.phi * A1.sA)
    G["X_b"] = -(A0.phi * A1b.sA) - A1.phi * A0b.sA \
        - A1b.h * A0.h - A0b.h * A1.h
    G["Wh_0"] = A0.sA * A1b.h - A0b.sA * A1.h
    G["Wh_1"] = A0.h * A1b.sA - A0b.h * A1.sA
    G["Wh_b"] = -(A0.h * A1b.sA) + A0b.sA * A1.h
    return G, AD, smooth


_GRP = ("Wphi", "Wp0", "Wp1", "X", "Wh")


def dred_grad_mv(B, T0, T1, branches=None):
    """Dred, its gradient (d/dth0, d/dth1, d/dbeta) and the D-atoms, with a
    mean-value refinement FIRST of each Wronskian group against its own
    partials, THEN of Dred against its assembled gradient.

    Returns (Dv, D0, D1, Db, AD, M, smooth, Dc)."""
    G, AD, smooth = _wc_groups(B, T0, T1, branches)
    bm = iv.mpf([mid(B), mid(B)])
    t0m = iv.mpf([mid(T0), mid(T0)])
    t1m = iv.mpf([mid(T1), mid(T1)])
    Gc, ADc, _ = _wc_groups(bm, t0m, t1m, branches)
    rb, r0, r1 = B - bm, T0 - t0m, T1 - t1m
    V = {}
    for nm in _GRP:
        V[nm] = _inter(G[nm], Gc[nm] + G[nm + "_0"] * r0 + G[nm + "_1"] * r1
                       + G[nm + "_b"] * rb)
    kap = PI_IV()
    HD, FD, sAD = AD.h, AD.phi, AD.sA
    M = kap * kap - FD * FD
    HDc, FDc = ADc.h, ADc.phi
    Mc = kap * kap - FDc * FDc
    # refine the D-atoms too (they are C^1 with hD' = sAD, phiD' = -hD;
    # D = th0 - th1 so d/dth0 = +, d/dth1 = -)
    HDr = _inter(HD, HDc + sAD * (r0 - r1))
    FDr = _inter(FD, FDc - HD * (r0 - r1))
    Mr = _inter(M, kap * kap - FDr * FDr)

    def assemble(hD, fD, m, W):
        return hD * hD * W["Wphi"] + hD * (kap * (W["Wp0"] + W["Wp1"])
                                           - fD * W["X"]) - m * W["Wh"]

    Dv = assemble(HDr, FDr, Mr, V)
    Dc = assemble(HDc, FDc, Mc, Gc)
    # gradient of Dred (using refined group values; sound: a valid enclosure
    # of each factor gives a valid enclosure of the derivative expression)
    Mp0 = 2 * FDr * HDr          # dM/dth0 ; dM/dth1 = -Mp0
    out = {}
    for var, sD in (("0", 1), ("1", -1)):
        dHD = sAD * sD
        dFD = -(HDr * sD)
        dM = Mp0 * sD
        out[var] = (2 * HDr * dHD * V["Wphi"] + HDr * HDr * G["Wphi_" + var]
                    + dHD * (kap * (V["Wp0"] + V["Wp1"]) - FDr * V["X"])
                    + HDr * (kap * (G["Wp0_" + var] + G["Wp1_" + var])
                             - dFD * V["X"] - FDr * G["X_" + var])
                    - dM * V["Wh"] - Mr * G["Wh_" + var])
    Db = (HDr * HDr * G["Wphi_b"]
          + HDr * (kap * (G["Wp0_b"] + G["Wp1_b"]) - FDr * G["X_b"])
          - Mr * G["Wh_b"])
    Dv = _inter(Dv, Dc + out["0"] * r0 + out["1"] * r1 + Db * rb)
    return Dv, out["0"], out["1"], Db, AD, Mr, smooth, Dc


def f_Dred(beta, t0, t1):
    """float mirror of Dred (diagnostic only)."""
    kap = PI
    D = t0 - t1
    FD, HD = f_phiJ(D), f_HJ(D)
    M = kap * kap - FD * FD
    f0, f1 = f_phiJ(t0), f_phiJ(t1)
    f0b, f1b = f_phiJ(t0 - beta), f_phiJ(t1 - beta)
    h0, h1 = f_HJ(t0), f_HJ(t1)
    h0b, h1b = f_HJ(t0 - beta), f_HJ(t1 - beta)
    Wphi = f1 * f0b - f0 * f1b
    Wp0 = f0 * h0b - f0b * h0
    Wp1 = f1 * h1b - f1b * h1
    X = f0 * h1b + f1 * h0b - f1b * h0 - f0b * h1
    Wh = h0 * h1b - h0b * h1
    return HD * HD * Wphi + HD * (kap * (Wp0 + Wp1) - FD * X) - M * Wh


def sep_det_grad_mv(B, T0, T1, branches=None):
    """Dsep, its gradient and (a0,a1), with the four column entries MEAN-VALUE
    refined against their own first partials.  Sound because a0,b0,a1,b1 are
    C^1 in (beta, th0, th1) (h and phi are C^1).  Returns
    (Dv, D0, D1, Db, a0, a1, AD, smooth, Dc) with Dc the thin-centre value."""
    C, AD, smooth = sep_columns(B, T0, T1, branches)
    bm = iv.mpf([mid(B), mid(B)])
    t0m = iv.mpf([mid(T0), mid(T0)])
    t1m = iv.mpf([mid(T1), mid(T1)])
    Cc, ADc, _ = sep_columns(bm, t0m, t1m, branches)
    rb, r0, r1 = B - bm, T0 - t0m, T1 - t1m
    V = {}
    for nm in ("a0", "b0", "a1", "b1"):
        V[nm] = _inter(C[nm], Cc[nm] + C[nm + "_0"] * r0 + C[nm + "_1"] * r1
                       + C[nm + "_b"] * rb)
    a0, b0, a1, b1 = V["a0"], V["b0"], V["a1"], V["b1"]
    Dv = a0 * b1 - a1 * b0
    D0 = C["a0_0"] * b1 + a0 * C["b1_0"] - C["a1_0"] * b0 - a1 * C["b0_0"]
    D1 = C["a0_1"] * b1 + a0 * C["b1_1"] - C["a1_1"] * b0 - a1 * C["b0_1"]
    Db = a0 * C["b1_b"] - C["a1_b"] * b0
    Dc = Cc["a0"] * Cc["b1"] - Cc["a1"] * Cc["b0"]
    Dv = _inter(Dv, Dc + D0 * r0 + D1 * r1 + Db * rb)
    return Dv, D0, D1, Db, a0, a1, AD, smooth, Dc
