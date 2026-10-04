"""
Centered two-student kernel: float mirror of website/widgets.js + rigorous
interval arithmetic (mpmath.iv) implementations.

Working precision is set by set_prec(); everything certified goes through
mpmath.iv (directed rounding).  float64 is used ONLY in the float mirror,
which exists to be compared bit-for-bit against widgets.js.
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
# float mirror of widgets.js (kernelOf("centered")) -- NOT part of the
# certified path.  Arithmetic order copied verbatim from the source lines.
# --------------------------------------------------------------------------


def f_mod(x, m):
    r = math.fmod(x, m)
    return r + m if r < 0 else r


def f_phiC(t):
    x = f_mod(t, PI)
    return (PI / 2 - x) * math.cos(x) + math.sin(x)


def f_HC(t):
    x = f_mod(t, PI)
    return (PI / 2 - x) * math.sin(x)


def f_dHC(t):
    x = f_mod(t, PI)
    return (PI / 2 - x) * math.cos(x) - math.sin(x)


def f_slopeAtom(t):
    return f_phiC(t) - 2 * abs(math.sin(t))


def f_Wtau(beta, t):
    return f_HC(t) * f_slopeAtom(t - beta) - f_HC(t - beta) * f_slopeAtom(t)


def f_Wpot(beta, t):
    return f_phiC(t) * f_HC(t - beta) - f_phiC(t - beta) * f_HC(t)


def f_Wwgt(beta, t):
    return abs(math.sin(t)) * f_HC(t - beta) - abs(math.sin(t - beta)) * f_HC(t)


# --------------------------------------------------------------------------
# interval helpers
# --------------------------------------------------------------------------

def I(a, b=None):
    if b is None:
        b = a
    return iv.mpf([a, b])


def PI_IV():
    return iv.pi


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


# --------------------------------------------------------------------------
# Branch-folded interval atoms.
#
# On the branch t in [k*pi, (k+1)*pi] put x = t - k*pi in [0, pi] and use the
# CLOSED BRANCH FORMS of widgets.js:
#     phi = (pi/2 - x) cos x + sin x        [phiC]
#     h   = (pi/2 - x) sin x                [HC   = sin t * arcsin(cos t)]
#     sA  = (pi/2 - x) cos x - sin x        [dHC  = phi - 2|sin| = h']
#     |sin t| = sin x
# and the branch derivatives (d/dt, dx/dt = 1)
#     phi' = -(pi/2 - x) sin x = -h
#     h'   = sA
#     sA'  = -2 cos x - (pi/2 - x) sin x
#     (|sin t|)' = cos x
# These are analytic on the OPEN branch; the kinks are exactly the branch
# endpoints x in {0, pi}.
# --------------------------------------------------------------------------

class Atoms(object):
    __slots__ = ("phi", "h", "sA", "asin", "dphi", "dh", "dsA", "dasin")


def atoms_on_branch(T, k):
    """T: interval of t known to lie in [k*pi, (k+1)*pi] (small overshoot ok)."""
    x = T - k * PI_IV()
    s = iv.sin(x)
    c = iv.cos(x)
    g = PI_IV() / 2 - x
    A = Atoms()
    A.phi = g * c + s
    A.h = g * s
    A.sA = g * c - s
    A.asin = s
    A.dphi = -(g * s)
    A.dh = A.sA
    A.dsA = -2 * c - g * s
    A.dasin = c
    return A


def endpoints(z):
    """Exact mp.mpf endpoints of an interval (mpmath's .a/.b are themselves
    zero-width intervals, which is useless for arithmetic)."""
    a = mpmath.mpf()
    a._mpf_ = z._mpi_[0]
    b = mpmath.mpf()
    b._mpf_ = z._mpi_[1]
    return a, b


def branch_pieces(T):
    """Rigorous cover of interval T by (piece, k) with piece inside branch k.

    Each piece is T intersected with an OUTWARD enclosure of [k*pi,(k+1)*pi];
    consecutive enclosures overlap, so the union of the pieces is all of T."""
    Ta, Tb = endpoints(T)
    klo = int(mpmath.floor(Ta / mpmath.mpf(3.14159265358979))) - 2
    khi = int(mpmath.floor(Tb / mpmath.mpf(3.14159265358980))) + 2
    out = []
    for k in range(klo, khi + 1):
        lo = endpoints(iv.mpf(k) * PI_IV())[0]          # <= k*pi
        hi = endpoints(iv.mpf(k + 1) * PI_IV())[1]      # >= (k+1)*pi
        plo = max(Ta, lo)
        phi_ = min(Tb, hi)
        if plo <= phi_:
            out.append((iv.mpf([plo, phi_]), k))
    assert out, "empty branch cover"
    # rigorous coverage check: the pieces must tile T with overlaps
    assert endpoints(out[0][0])[0] == Ta and endpoints(out[-1][0])[1] == Tb
    for i in range(len(out) - 1):
        assert endpoints(out[i + 1][0])[0] <= endpoints(out[i][0])[1]
    return out


def atoms(T):
    """Enclosure of all atoms over T, valid across kinks (hull over branches)."""
    ps = branch_pieces(T)
    As = [atoms_on_branch(p, k) for p, k in ps]
    R = Atoms()
    for f in Atoms.__slots__:
        R.__setattr__(f, hull(*[getattr(A, f) for A in As]))
    return R, (len(ps) == 1)


# --------------------------------------------------------------------------
# The three mass-free determinants, mirrored from the widgets.js atoms.
# Derivative in t uses the exact cancellations
#     d/dt Wtau = h(t) sA'(t-b) - h(t-b) sA'(t)
#     d/dt Wpot = phi(t) sA(t-b) - phi(t-b) sA(t)
# (the other two terms cancel identically), valid on a smooth branch pair.
# --------------------------------------------------------------------------

def dets_branch(B, T, kt, ku):
    """Same as dets_iv but with the branch indices supplied by the caller.

    Legitimate exactly when T is contained in the CLOSED branch [kt*pi,(kt+1)*pi]
    and T-B in [ku*pi,(ku+1)*pi]: on a closed branch the closed forms agree with
    the periodic functions (the kink sits at the branch endpoint, where the
    one-sided closed form still gives the correct value).  This is how the kink
    lattice is handled: it becomes the boundary of the smooth pieces."""
    At = atoms_on_branch(T, kt)
    Au = atoms_on_branch(T - B, ku)
    Wt = At.h * Au.sA - Au.h * At.sA
    Wp = At.phi * Au.h - Au.phi * At.h
    Ww = At.asin * Au.h - Au.asin * At.h
    dWt = At.h * Au.dsA - Au.h * At.dsA
    dWp = At.phi * Au.sA - Au.phi * At.sA
    dWw = At.dasin * Au.h + At.asin * Au.dh - Au.dasin * At.h - Au.asin * At.dh
    return {"Wtau": (Wt, dWt), "Wpot": (Wp, dWp), "Wwgt": (Ww, dWw)}


def dets_iv(B, T):
    """B, T intervals.  Returns dict name -> (F, dF/dt, smooth_flag)."""
    At, st_ok = atoms(T)
    Au, su_ok = atoms(T - B)
    sm = st_ok and su_ok
    Wt = At.h * Au.sA - Au.h * At.sA
    Wp = At.phi * Au.h - Au.phi * At.h
    Ww = At.asin * Au.h - Au.asin * At.h
    dWt = At.h * Au.dsA - Au.h * At.dsA
    dWp = At.phi * Au.sA - Au.phi * At.sA
    dWw = At.dasin * Au.h + At.asin * Au.dh - Au.dasin * At.h - Au.asin * At.dh
    return {"Wtau": (Wt, dWt, sm), "Wpot": (Wp, dWp, sm), "Wwgt": (Ww, dWw, sm)}


# --------------------------------------------------------------------------
# REDUCED (kink-free) piecewise forms.
#
# For beta in (0, pi) the period [0, pi] of t splits at the kink lattice into
#   piece A : t in (beta, pi)   ->  b = beta,      sigma = +1, w = t - (beta+pi)/2
#   piece B : t in (0, beta)    ->  b = pi - beta, sigma = -1, w = t - beta/2
# and on each piece, with W = (pi - b)/2 and |w| <= W,
#   Wtau = sigma * Phi(b,w),  Wpot = sigma * Psi(b,w),  Wwgt = sigma * Om(b,w)
#   Q  (b,w) = (w^2 - b^2/4) sin b
#   Om (b,w) = (b/2)(cos 2w + cos b)  = b cos(w + b/2) cos(w - b/2)
#   Phi = Q + Om,    Psi = -Q + Om      [hence Phi + Psi = 2 Om identically]
# Both Phi and Psi are ANALYTIC in (b,w) on the whole piece: the reduction has
# absorbed the kinks into the piece boundaries |w| = W.
# --------------------------------------------------------------------------

def redu(b, w):
    """b, w intervals -> (Phi, Psi, Om, Phi_w, Psi_w, Om_w)."""
    sb = iv.sin(b)
    cb = iv.cos(b)
    Q = (w * w - b * b / 4) * sb
    Om = (b / 2) * (iv.cos(2 * w) + cb)
    s2w = iv.sin(2 * w)
    Qw = 2 * w * sb
    Omw = -b * s2w
    return (Q + Om, -Q + Om, Om, Qw + Omw, -Qw + Omw, Omw)


def piece_params(beta_iv, piece):
    """piece in {'A','B'} -> (b, sigma)."""
    if piece == "A":
        return beta_iv, 1
    return PI_IV() - beta_iv, -1


def t_from_w(beta_iv, piece, w):
    if piece == "A":
        return (beta_iv + PI_IV()) / 2 + w
    return beta_iv / 2 + w


# --------------------------------------------------------------------------
# mass direction / y coordinate of a torque root
#   (s0, s1) proportional to (-h(t-beta), h(t)),  psi = atan2(s0,s1) mod pi,
#   y = 2 psi / pi - 1  in (-1, 1]              [widgets.js ratioCoord]
# --------------------------------------------------------------------------

def y_of_root(B, T):
    At, _ = atoms(T)
    Au, _ = atoms(T - B)
    s0 = -Au.h
    s1 = At.h
    return s0, s1


def y_from_s(s0, s1):
    """Certified y enclosure from an (s0,s1) enclosure.

    y = 2*psi/pi - 1 with psi = atan2(s0,s1) folded by the (-,-) isometry into
    (0, pi] -- exactly widgets.js ratioCoord.  Returns None if the enclosure
    straddles the fold seam (then y is not determined by this box)."""
    if contains_zero(s0) and contains_zero(s1):
        return None
    a = iv.atan2(s0, s1)                      # subset of (-pi, pi]
    if lo(a) > 0:
        psi = a
    elif hi(a) < 0:
        psi = a + PI_IV()
    else:
        return None
    return 2 * psi / PI_IV() - 1
