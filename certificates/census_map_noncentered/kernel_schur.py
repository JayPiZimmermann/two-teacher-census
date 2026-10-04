"""The MASS-FREE separated type-boundary function, in interval arithmetic.

This is the one remaining input of the face-level count statement.  The Lean
side (`CountConstancyJ/FoldFreeInputJ.censusAngleMapJ_foldFree_of_kernelTeacher
SchurDetT`) says: at a zero of the census angle map, off the coincidence and
opposition lattices and with both teacher torques alive, nonvanishing of
`generalJKernelTeacherSchurDetT` gives the fold-freeness the counting schema
consumes.  It carries NO teacher mass -- its arguments are the teacher gap and
the two student angles -- which is exactly why one certificate can serve a
whole face of the census, whose `y` coordinate only moves the masses.

TRANSCRIBED FROM THE TREE, definition by definition, so a reader can check it
against the Lean rather than against a paraphrase
(`NoncenteredGeneral/SchurBoundary/SeparatedTypeBoundaryJ.lean`,
`SeparatedCount/ScalarColumnJ.lean`):

    num(D,x)  = (pi*phi(x-D) - phi(D)*phi(x))*h(D) - (pi^2-phi(D)^2)*h(x)
    S0        = -num(D, th0 - beta)          [ = -generalJAngleEq0 beta 0 1 ]
    S1        =  num(D, th0)                 [ =  generalJAngleEq0 beta 1 0 ]

the KERNEL VECTOR of the first eliminated angle equation, and then the Cramer
Schur entries at that vector, with `V(g) = S0*phi(g) + S1*phi(g-beta)`,
`L(g) = S0*|sin g| + S1*|sin(g-beta)|`, `M = pi^2 - phi(D)^2`,

    N0 = pi*V(th0) - phi(D)*V(th1),   N1 = pi*V(th1) - phi(D)*V(th0)
    K  = 2*|sin D| - phi(D)
    T00 = M*N0*N1*K - M^2*N0*(2*L(th0) - V(th0)) - pi*h(D)^2*N0^2
    T11 = M*N0*N1*K - M^2*N1*(2*L(th1) - V(th1)) - pi*h(D)^2*N1^2
    T01 = -(M*N0*N1*K) - phi(D)*h(D)^2*N0*N1
    detT = T00*T11 - T01^2

`|sin t|` is taken as `sin` of the branch coordinate, which lies in `[0, pi]`,
so it is exact rather than an absolute value of an interval straddling zero.
"""
import mpmath

import census_cert as X
import noncentered as J

PIv = X.PIv


def _abs_sin(T):
    """|sin| over an interval, with the half-branch orientation included.

    On branch ``n``, ``p`` is nonnegative for even ``n`` and nonpositive for
    odd ``n``.  Hence ``|sin T| = e*sin(p)``, not always ``sin(p)``.  The old
    implementation omitted ``e`` and returned the negative of the required
    atom on every odd branch.
    """
    ps = J.branch_pieces(T)
    out = None
    for piece, n in ps:
        e = 1 if n % 2 == 0 else -1
        v = e * X.iv.sin(J.p_of(piece, n))
        out = v if out is None else J.hull(out, v)
    return out


def _atoms(T):
    A, _ = J.atoms(T)
    return A


def num(D, x):
    """`separatedNumJ D x` of the tree, in interval arithmetic."""
    AD = _atoms(D)
    Ax = _atoms(x)
    Axd = _atoms(x - D)
    M = PIv * PIv - AD.phi * AD.phi
    return (PIv * Axd.phi - AD.phi * Ax.phi) * AD.h - M * Ax.h


def kernel_teacher_schur_detT(beta, th0, th1):
    """`generalJKernelTeacherSchurDetT beta ![th0, th1]`, interval."""
    D = th0 - th1
    AD = _atoms(D)
    phiD, hD = AD.phi, AD.h
    M = PIv * PIv - phiD * phiD
    S0 = -num(D, th0 - beta)
    S1 = num(D, th0)

    def V(g):
        return S0 * _atoms(g).phi + S1 * _atoms(g - beta).phi

    def L(g):
        return S0 * _abs_sin(g) + S1 * _abs_sin(g - beta)

    V0, V1 = V(th0), V(th1)
    N0 = PIv * V0 - phiD * V1
    N1 = PIv * V1 - phiD * V0
    K = 2 * _abs_sin(D) - phiD
    base = M * N0 * N1 * K
    T00 = base - M * M * N0 * (2 * L(th0) - V0) - PIv * hD * hD * N0 * N0
    T11 = base - M * M * N1 * (2 * L(th1) - V1) - PIv * hD * hD * N1 * N1
    T01 = -base - phiD * hD * hD * N0 * N1
    return T00 * T11 - T01 * T01, (T00, T11, T01, M, N0, N1)
