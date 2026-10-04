"""REGRESSION GUARD for the separated residual's interval Jacobian.

WHY THIS FILE EXISTS (2026-08-08).  `census_cert.sep_res(..., want_jac=True)`
returned a Jacobian block whose `dc1_1` entry carried `+ hD * P0` where the
chain rule gives `- hD * P0`.  The consequence was NOT a wide enclosure but
an INVALID one: the interval `j01` did not contain the true range of
`dG0/d(th1)` in 188 of 300 random boxes.  An invalid derivative enclosure
breaks the Krawczyk test's hypothesis, and the test then returned spurious
`unique` verdicts -- i.e. it CERTIFIED SEPARATED FAMILIES THAT DO NOT EXIST.
At the F4 witness (beta, y) = (0.30, -0.50) that fabricated 8 mutually
disjoint "families" in a strip where a 241x241 float scan puts
min(|G0| + |G1|) at 1.5e-2 and Newton escapes the box.

The bug was invisible to every check the campaign had:
  * the RESIDUAL is correct -- `census_cert.sep_res` and the `count_scan`
    float mirror agree to 9 digits, so residual-level validation passes;
  * the interval layer is sound (mpmath.iv rounds outward), so precision
    sweeps say nothing -- a WIDER interval makes Krawczyk FAIL, never
    succeed spuriously.  Only an enclosure that is WRONG can manufacture a
    verdict, and only a derivative check can see it;
  * `n_separated` counts looked plausible and matched a census that had
    been produced by the same machinery.
The check that DID catch it is the one below, and the shape generalises:
**an interval enclosure of a DERIVATIVE must be validated against sampled
finite differences of the function it claims to bound, not against another
evaluation of itself.**  (Formalization skill 1j: two implementations sharing
a model verify consistency, not soundness.)

The audit: over random teachers and random boxes spanning the enumeration
domain, sample the true partials by finite differences of the independent
float mirror and assert each is inside the interval Jacobian.

Validity domain: float64 finite differences at h = 1e-7 with a 1e-5 slack,
boxes of side 0.02 kept off the branch lattice, teachers with both masses
bounded away from zero.  The audit can only FALSIFY an enclosure; it cannot
certify one.  Exit 0 = no violation found.

Usage:  python3 sep_jacobian_audit.py [n_boxes]
"""
import math
import sys

import mpmath
import numpy as np

import census_cert as X
import count_scan as CS
import noncentered as J

H = 1e-7
SLACK = 1e-5
SIDE = 0.02
NAMES = ("j00 = dG0/dth0", "j01 = dG0/dth1",
         "j10 = dG1/dth0", "j11 = dG1/dth1")


def audit(n_boxes=300, seed=1, verbose=False):
    rng = np.random.default_rng(seed)
    bad = [0, 0, 0, 0]
    worst = [0.0, 0.0, 0.0, 0.0]
    tested = 0
    examples = []
    while tested < n_boxes:
        bf = round(float(rng.uniform(0.2, 3.0)), 4)
        yf = round(float(rng.uniform(-0.9, 0.9)), 4)
        B, Y = J.I(repr(bf)), J.I(repr(yf))
        s0i, s1i = X.masses_at(Y)
        Tc = X.Teacher(B, s0i, s1i)
        psi = (yf + 1) * math.pi / 2
        s0, s1 = math.sin(psi), math.cos(psi)

        sa = float(rng.uniform(0.05, 6.0))
        da = float(rng.uniform(0.10, 2.90))
        box = (sa, sa + SIDE, da, da + SIDE)
        T1 = X.iv.mpf([mpmath.mpf(repr(box[0])), mpmath.mpf(repr(box[1]))])
        Dm = X.iv.mpf([mpmath.mpf(repr(box[2])), mpmath.mpf(repr(box[3]))])

        gs = np.linspace(box[0], box[1], 25)
        gd = np.linspace(box[2], box[3], 25)
        S, D = np.meshgrid(gs, gd, indexing="ij")
        S, D = S.ravel(), D.ravel()
        f0, f1 = CS.v_res(bf, s0, s1, S + D, S)
        a0, a1 = CS.v_res(bf, s0, s1, S + D + H, S)
        b0, b1 = CS.v_res(bf, s0, s1, S + D, S + H)
        samp = [(a0 - f0) / H, (b0 - f0) / H, (a1 - f1) / H, (b1 - f1) / H]
        if not all(np.all(np.isfinite(p)) for p in samp):
            continue

        r = X.sep_res(Tc, T1 + Dm, T1, None, want_jac=True, Dm=Dm)
        if r is None:
            continue
        tested += 1
        for k, q in enumerate(r[8]):
            lo, hi = float(J.lo(q)), float(J.hi(q))
            p = samp[k]
            miss = max(lo - float(np.nanmin(p)), float(np.nanmax(p)) - hi)
            if miss > SLACK:
                bad[k] += 1
                worst[k] = max(worst[k], miss)
                if len(examples) < 5:
                    examples.append(
                        {"entry": NAMES[k], "beta": bf, "y": yf,
                         "box": box, "enclosure": [lo, hi],
                         "sampled": [float(np.nanmin(p)),
                                     float(np.nanmax(p))],
                         "miss": miss})
        if verbose and tested % 50 == 0:
            print("  %d boxes, violations %s" % (tested, bad), flush=True)
    return tested, bad, worst, examples


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    tested, bad, worst, examples = audit(n, verbose=True)
    print("boxes tested:", tested)
    for k in range(4):
        print("  %-16s invalid enclosures: %4d   worst miss: %.3e"
              % (NAMES[k], bad[k], worst[k]))
    for e in examples:
        print("  EXAMPLE", e)
    if any(bad):
        print("SEP JACOBIAN AUDIT FAILED -- the Krawczyk test's hypothesis "
              "is violated, so every 'unique'/'empty' verdict downstream is "
              "void.")
        sys.exit(1)
    print("SEP JACOBIAN AUDIT PASSED (%d boxes, no invalid enclosure)" % tested)


if __name__ == "__main__":
    main()
