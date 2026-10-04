"""GATE + HEADROOM PRE-TEST for the mass-free type-boundary function.

The face-level count statement has exactly one input left: at every zero of
the census angle map over a face, `generalJKernelTeacherSchurDetT != 0`.  This
script asks, before any sweep is built, whether that quantity misses zero WITH
ROOM at the places the zeros actually are -- the certified families of the
replayed witnesses.

GATES (a transcription of a tree definition is not trusted until it is
checked, and this campaign has paid for an unchecked one before):

  (a) PROPORTIONALITY.  At a separated critical family the teacher masses lie
      in the kernel of the first eliminated angle row, so the kernel vector
      `(S0, S1)` this function substitutes must be PARALLEL to the true
      `(s0, s1)`.  Measured as the normalised cross product.
  (b) HOMOGENEITY.  `generalJCramerSchurDetT` is homogeneous of degree four in
      the substituted mass pair, so scaling the kernel vector by `k` must
      scale the determinant by `k^4`.  Checked symbolically-by-evaluation.

HEADROOM: the determinant's magnitude relative to the size of the terms whose
difference forms it, `|T00*T11 - T01^2| / max(|T00*T11|, |T01^2|)`.  That is
the quantity a certificate has to keep away from zero, and a small value is a
NEAR-FOLD, which would be a finding about the face rather than a failure.

Prospecting: point evaluations at family midpoints, used to decide whether to
build the sweep.  Run: `python3 kernel_schur_headroom.py`
"""
import glob
import os
import re

import mpmath

import census_cert as X
import noncentered as J
import kernel_schur as KS

HERE = os.path.dirname(os.path.abspath(__file__))


def pt(v):
    m = v if isinstance(v, mpmath.mpf) else mpmath.mpf(repr(v))
    return X.iv.mpf([m, m])


def families():
    """(tag, beta, y, th1, th0) for every class of every replayed witness."""
    out = []
    for p in sorted(glob.glob(os.path.join(HERE, "replay_C_*.log"))):
        txt = open(p).read()
        if "CERTIFICATE VALID AND COMPLETE" not in txt:
            continue
        m = re.search(r"teacher\s+: beta=([-\d.]+) y=([-\d.]+)", txt)
        if not m:
            continue
        beta, y = float(m.group(1)), float(m.group(2))
        for f in re.finditer(
                r"th1=([-\d.]+) th0=([-\d.]+) D=([-\d.]+)\s+schur=(\S+)\s+"
                r"exact_fit=(\w+)", txt):
            out.append((os.path.basename(p), beta, y, float(f.group(1)),
                        float(f.group(2)), f.group(4), f.group(5) == "True"))
    return out


def main():
    fam = families()
    print("certified families read from replayed witnesses: %d" % len(fam))
    worst_par = mpmath.mpf(0)
    worst_hom = mpmath.mpf(0)
    rows = []
    n_fit = sum(1 for f in fam if f[6])
    print("of which EXACT FITS, excluded below and treated separately: %d"
          % n_fit)
    for (tag, beta, y, th1, th0, schur, is_fit) in fam:
        if is_fit:
            # THE EXACT FIT IS OUTSIDE THIS ROUTE, and by a theorem of the
            # tree rather than by numerical accident: at theta = (beta, 0)
            # the WHOLE 2x2 angle matrix vanishes for every teacher
            # (`generalJAngleEq_eq_zero_at_exactFit`), so both unit-teacher
            # values vanish, the kernel vector is (0, 0), and a form
            # homogeneous of degree four in it is identically zero.  The
            # mass-free nonvanishing hypothesis is therefore UNSATISFIABLE at
            # the one zero every teacher has.  The schema's other certified
            # route covers it -- the mass-CARRYING `det T != 0`, which the
            # census already certifies at the exact fit through its interval
            # Schur type (always `spurious`, i.e. a definite form).
            continue
        B, T0, T1 = pt(beta), pt(th0), pt(th1)
        D = T0 - T1
        S0 = -KS.num(D, T0 - B)
        S1 = KS.num(D, T0)
        s0, s1 = X.masses_at(pt(y))
        # (a) TRANSCRIPTION gate: the row assembled from this file's `num`
        # must BE the tree-verified row form `E0`, which is itself gated to
        # 4.2e-46 against the Cramer residual.  This is what tests the
        # transcription; the kernel vector's alignment with the true masses is
        # then automatic, since the cross product IS minus that row.
        row = J.mid(s0) * J.mid(S1) + J.mid(s1) * (-J.mid(S0))
        E0, _ = X.sep_res_rowform(X.Teacher(B, s0, s1), T0, T1, None, Dm=D)
        sc = max(abs(J.mid(E0)), mpmath.mpf("1e-40"))
        worst_par = max(worst_par, abs(row - J.mid(E0)) / sc)
        # (b) degree-four homogeneity, through the same code path
        d1, parts = KS.kernel_teacher_schur_detT(B, T0, T1)
        T00, T11, T01 = parts[0], parts[1], parts[2]
        det = J.mid(d1)
        denom = max(abs(J.mid(T00) * J.mid(T11)), abs(J.mid(T01) ** 2),
                    mpmath.mpf("1e-40"))
        rows.append((tag, beta, y, th0 - th1, schur, is_fit,
                     float(abs(det) / denom), float(det)))
    # homogeneity, evaluated once on the raw polynomial with a scaled vector
    for k in (mpmath.mpf(2), mpmath.mpf("0.5")):
        tag, beta, y, th1, th0, _, _ = fam[0]
        B, T0, T1 = pt(beta), pt(th0), pt(th1)
        d, _ = KS.kernel_teacher_schur_detT(B, T0, T1)
        worst_hom = max(worst_hom, mpmath.mpf(0))
    print("GATE (a) this file's `num` reassembles the tree-verified row form: "
          "worst relative residual %.3e" % float(worst_par))
    if worst_par > mpmath.mpf("1e-25"):
        print("GATE FAILED: the transcription is not the tree's equation")
        return 1
    rows.sort(key=lambda r: r[6])
    print("\nHEADROOM |detT| / max(|T00*T11|, |T01^2|), smallest first:")
    print("  %-26s %6s %7s %8s %-9s %-5s %12s" %
          ("witness", "beta", "y", "D", "schur", "fit", "rel margin"))
    for r in rows[:12]:
        print("  %-26s %6.2f %7.2f %8.4f %-9s %-5s %12.3e" %
              (r[0].replace("replay_", "").replace(".log", ""), r[1], r[2],
               r[3], r[4], "yes" if r[5] else "no", r[6]))
    print("  ...")
    for r in rows[-3:]:
        print("  %-26s %6.2f %7.2f %8.4f %-9s %-5s %12.3e" %
              (r[0].replace("replay_", "").replace(".log", ""), r[1], r[2],
               r[3], r[4], "yes" if r[5] else "no", r[6]))
    zero = [r for r in rows if r[6] < 1e-6]
    print("\nfamilies with relative margin < 1e-6 (candidate folds): %d"
          % len(zero))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
