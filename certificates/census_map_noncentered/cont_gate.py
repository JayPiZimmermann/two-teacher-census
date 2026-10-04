"""GATE for the angle-map Jacobian determinant used by `cont_cert.py`.

A FALSIFIER, not evidence -- exactly like `validate_rowform.py` and
`sep_jacobian_audit.py`.  It cannot certify the identity; it can refute it, and
an invalid derivative enclosure is a failure this certificate has already paid
for once (the 2026-08-08 Krawczyk sign fix), so the quantity is gated before
any leaf is allowed to claim fold-freeness with it.

THREE CHECKS.

(1) THE CHART CHANGE.  `cont_cert.jac_det_theta` computes

        det_theta = dE0/dD * dE1/ds - dE0/ds * dE1/dD

from `census_cert.sep_res_rowform_jac`, whose derivatives are in the `(s, D)`
chart `th0 = s + D`, `th1 = s`.  Since `d/dth0 = d/dD` and
`d/dth1 = d/ds - d/dD`, that expression is the determinant of the Jacobian of
`(E0, E1)` in the `(th0, th1)` chart.  Checked against central differences of
`census_cert.sep_res_rowform` in `th0` and `th1` at point boxes.

(2) THE TREE'S OWN CLOSED FORMS.  `LeanFormalization/Planar/SignedN2/FreeMassJ/
SeparatedIndex/AngleJacJ.lean` DEFINES the four partial derivatives as closed
forms and PROVES them to be the derivatives of `generalJAngleEq0/1`.  Those
four are transcribed here verbatim and their determinant
`separatedAngleJacDetJ` is compared with `-det_theta`.  The sign is not a
convention to be chosen: the tree's second equation is

    generalJAngleEq1 = generalJMassNum0 * couplingHJ (th1 - th0)
                       + generalJMassDet * generalTorqueJ th1

and `couplingHJ` is ODD, so it is the NEGATIVE of this file's `E1`
(`E1 = (pi*P0 - phi(D)*P1)*h(D) - massDet*A1`).  One row negated flips the
determinant, hence the minus.  Both are evaluated at the SAME exact `(th0,
th1)`, since rounding `th0 = s + D` in floating point alone shows up at
`1e-16` and would mask a real defect.

(3) THE TRANSCRIPTION AGAINST ITS OWN EQUATIONS.  The four transcribed closed
forms are checked against central differences of the transcribed
`generalJAngleEq0/1`, so a typo in (2) cannot be absorbed by a matching typo
in the determinant.

Run: `python3 cont_gate.py`
"""
import random
import sys

import mpmath
from mpmath import iv

import census_cert as X
import cont_cert as C
import noncentered as J

random.seed(20260809)


def point(v):
    m = v if isinstance(v, mpmath.mpf) else mpmath.mpf(repr(v))
    return iv.mpf([m, m])


# --------------------------------------------------------------------------
# transcription of AngleJacJ.lean
# --------------------------------------------------------------------------

def lean_entries(Tc, t0, t1):
    """`separatedAngleJac00J/01J/10J/11J`, transcribed.  Vocabulary:
    `generalTeacherPotentialJ = P`, `generalTorqueJ = A`,
    `generalJTorqueSlope = A'`, `phiCosJ = K`, `couplingHJ = h`, and the
    coupling slope `K(D) - 2|sin D|`, which is the tree's `sA`."""
    T0, T1 = point(t0), point(t1)
    Dm = T0 - T1
    AD, _ = J.atoms(Dm)
    K, h, hp = AD.phi, AD.h, AD.sA
    kap = J.PI_IV()
    P0, A0, Ad0 = Tc.load(T0)
    P1, A1, Ad1 = Tc.load(T1)
    massDet = kap * kap - K * K
    num0 = kap * P0 - K * P1
    num1 = kap * P1 - K * P0
    j00 = (h * P0 - K * A0) * h + num1 * hp + 2 * K * h * A0 + massDet * Ad0
    j01 = (kap * A1 - h * P0) * h - num1 * hp - 2 * K * h * A0
    j10 = -((kap * A0 + h * P1) * h) - num0 * hp + 2 * K * h * A1
    j11 = (h * P1 + K * A1) * h + num0 * hp - 2 * K * h * A1 + massDet * Ad1
    return j00, j01, j10, j11


def lean_eqs(Tc, t0, t1):
    """`generalJAngleEq0`, `generalJAngleEq1`, transcribed."""
    T0, T1 = point(t0), point(t1)
    Dm = T0 - T1
    AD, _ = J.atoms(Dm)
    K, h = AD.phi, AD.h
    kap = J.PI_IV()
    P0, A0, _ = Tc.load(T0)
    P1, A1, _ = Tc.load(T1)
    massDet = kap * kap - K * K
    eq0 = (kap * P1 - K * P0) * h + massDet * A0
    eq1 = -((kap * P0 - K * P1) * h) + massDet * A1
    return eq0, eq1


def rand_teacher():
    b = random.uniform(0.15, 3.0)
    y = random.uniform(-0.97, -0.03)
    B = J.I(str(b))
    Y = J.I(str(y))
    s0, s1 = X.masses_at(Y)
    return X.Teacher(B, s0, s1), b, y


def main():
    J.set_prec(300)
    n = 200
    w1 = w2 = w3 = mpmath.mpf(0)
    hstep = mpmath.mpf(10) ** (-25)
    for _ in range(n):
        Tc, b, y = rand_teacher()
        s = mpmath.mpf(repr(random.uniform(0.05, 6.2)))
        D = mpmath.mpf(repr(random.uniform(0.25, 3.05)))
        t0, t1 = s + D, s                      # EXACT, no float rounding
        S, DD = point(s), point(D)

        det_theta = C.jac_det_theta(Tc, S + DD, S, None, DD)

        # (1) against central differences of the row form in (th0, th1)
        def E(a, bq):
            return X.sep_res_rowform(Tc, point(a), point(bq), None,
                                     Dm=point(a) - point(bq))
        pa, pb = E(t0 + hstep, t1), E(t0 - hstep, t1)
        pc, pd = E(t0, t1 + hstep), E(t0, t1 - hstep)
        f00 = J.mid((pa[0] - pb[0]) / (2 * hstep))
        f01 = J.mid((pc[0] - pd[0]) / (2 * hstep))
        f10 = J.mid((pa[1] - pb[1]) / (2 * hstep))
        f11 = J.mid((pc[1] - pd[1]) / (2 * hstep))
        fd = f00 * f11 - f01 * f10
        w1 = max(w1, abs(J.mid(det_theta) - fd) / max(abs(fd), mpmath.mpf(1)))

        # (2) against the tree's closed forms
        j00, j01, j10, j11 = lean_entries(Tc, t0, t1)
        lean_det = J.mid(j00 * j11 - j01 * j10)
        w2 = max(w2, abs(lean_det + J.mid(det_theta))
                 / max(abs(lean_det), mpmath.mpf(1)))

        # (3) the transcription against its own equations
        qa, qb = lean_eqs(Tc, t0 + hstep, t1), lean_eqs(Tc, t0 - hstep, t1)
        qc, qd = lean_eqs(Tc, t0, t1 + hstep), lean_eqs(Tc, t0, t1 - hstep)
        g = [J.mid((qa[0] - qb[0]) / (2 * hstep)),
             J.mid((qc[0] - qd[0]) / (2 * hstep)),
             J.mid((qa[1] - qb[1]) / (2 * hstep)),
             J.mid((qc[1] - qd[1]) / (2 * hstep))]
        for k, jj in enumerate((j00, j01, j10, j11)):
            w3 = max(w3, abs(J.mid(jj) - g[k]) / max(abs(g[k]), mpmath.mpf(1)))

    print("samples                                        : %d" % n)
    print("(1) det_theta vs central differences of E       : %s"
          % mpmath.nstr(w1, 4))
    print("(2) det_theta + separatedAngleJacDetJ (tree)    : %s"
          % mpmath.nstr(w2, 4))
    print("(3) tree closed forms vs their own equations    : %s"
          % mpmath.nstr(w3, 4))
    tolfd = mpmath.mpf(10) ** (-18)
    tolex = mpmath.mpf(10) ** (-40)
    bad = (w1 > tolfd) or (w3 > tolfd) or (w2 > tolex)
    print("GATE %s (finite-difference tolerance 1e-18, exact-identity "
          "tolerance 1e-40)" % ("FAILED" if bad else "PASSED"))
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
