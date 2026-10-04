"""GATE for the row-7a certificate's MAP verdict.

The MAP verdict kills a box when an enclosure of `census_cert.sep_res_rowform`
misses zero, and it is sound only if that pair really is the census angle map's
pair of eliminated angle equations — otherwise the certificate would be
excluding zeros of a different function.  Two links are already gated
elsewhere and one was not:

  * `validate_rowform.py` — `E_i = massDet(D)·G_i` exactly (4.2e-46), so a
    nonzero `E_i` gives a nonzero `G_i`, the direction the exclusion needs;
  * `kernel_schur_headroom.py` gate (a) — the FIRST row assembled from the
    transcribed `separatedNumJ` reproduces `E0`;
  * THIS FILE — the SECOND row too.  `cont_gate.py` tests the row form's
    Jacobian determinant against `separatedAngleJacDetJ`, which would very
    probably catch a wrong `E1`, but "very probably" is not a gate.

Both rows are assembled from `kernel_schur.num` (the transcription of the
tree's `separatedNumJ`) at the arguments the tree's own unit-teacher lemmas
use, and compared UP TO SIGN, since a global sign cannot change whether an
enclosure misses zero.

Run: `python3 row7a_rowgate.py`
"""
import random

import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS

random.seed(20260810)


def pt(v):
    m = mpmath.mpf(repr(v))
    return X.iv.mpf([m, m])


def main():
    worst = [mpmath.mpf(0), mpmath.mpf(0)]
    n = 0
    for _ in range(200):
        b = random.uniform(0.2, 3.0)
        s = random.uniform(0.3, 6.0)
        d = random.uniform(0.1, 3.0)
        y = random.uniform(-0.95, 0.95)
        B, T1, Dm = pt(b), pt(s), pt(d)
        T0 = T1 + Dm
        s0, s1 = X.masses_at(pt(y))
        Tc = X.Teacher(B, s0, s1)
        E0, E1 = X.sep_res_rowform(Tc, T0, T1, None, Dm=Dm)
        # student 0's two teacher offsets; student 1's two teacher offsets
        r0 = (J.mid(s0) * J.mid(KS.num(Dm, T0))
              + J.mid(s1) * J.mid(KS.num(Dm, T0 - B)))
        r1 = (J.mid(s0) * J.mid(KS.num(Dm, -T1))
              + J.mid(s1) * J.mid(KS.num(Dm, B - T1)))
        for i, (r, E) in enumerate(((r0, E0), (r1, E1))):
            sc = max(abs(J.mid(E)), mpmath.mpf("1e-30"))
            worst[i] = max(worst[i],
                           min(abs(r - J.mid(E)), abs(r + J.mid(E))) / sc)
        n += 1
    print("row 0 vs sep_res_rowform E0 (up to sign): %d samples, worst "
          "relative residual %.3e" % (n, float(worst[0])))
    print("row 1 vs sep_res_rowform E1 (up to sign): worst relative residual "
          "%.3e" % float(worst[1]))
    if max(worst) > mpmath.mpf("1e-25"):
        print("GATE FAILED: the MAP verdict is not testing the census angle "
              "map")
        return 1
    print("gate passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
