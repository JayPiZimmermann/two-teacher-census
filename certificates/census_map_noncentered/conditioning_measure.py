"""Measured conditioning of the two-dimensional mass-free exclusion sweep
(the artifact behind the sweep's honest verdict in CERTIFICATE.md).

At a THIN beta the half domain (s, D) = (th1, th0-th1), s in [0, 2pi],
D in [delta, pi - delta], is covered by uniform boxes of side w and each box is
tested by the exclusion Psi = Dsep / h(D)^4 (mean-value refined, simplified
gradient).  Reported per w: the number of boxes and the fraction excluded.

At the finest w, additionally:
  * the interval-overestimation ratios  width(enclosure(Dsep)) / true range
    and the same for dDsep/dth0, on a subsample of boxes (true range estimated
    by a dense float evaluation on the box);
  * the fraction of cells whose float sign map shows a genuine sign change
    (an estimate of how many cells actually contain the curve).

Writes conditioning.json.
"""
import json
import math
import os
import time

import mpmath
from mpmath import iv

import noncentered as J
import certify_separated as S

HERE = os.path.dirname(os.path.abspath(__file__))
J.set_prec(S.PREC)

BETA = mpmath.mpf("1.22")
DELTA = mpmath.mpf("0.15")
B = iv.mpf([BETA, BETA])
two_pi = float(J.mid(2 * J.PI_IV()))
dpi = float(J.mid(J.PI_IV()))

out = {"beta": float(BETA), "delta": float(DELTA),
       "precision_bits": S.PREC, "rows": []}


def boxes_of(w):
    ns = int(math.ceil(two_pi / w))
    d_lo, d_hi = float(DELTA), dpi - float(DELTA)
    nd = int(math.ceil((d_hi - d_lo) / w))
    bx = []
    for i in range(ns):
        for k in range(nd):
            bx.append((two_pi * i / ns, two_pi * (i + 1) / ns,
                       d_lo + (d_hi - d_lo) * k / nd,
                       d_lo + (d_hi - d_lo) * (k + 1) / nd))
    return bx


t0 = time.time()
finest = None
for w in (0.1, 0.05, 0.025):
    bx = boxes_of(w)
    nexc = 0
    for (s0, s1, d0, d1) in bx:
        T1 = iv.mpf([s0, s1])
        Dm = iv.mpf([d0, d1])
        T0 = T1 + Dm
        br = S.branches_for(B, T0, T1)
        r = S.evaluate(B, T0, T1, br)
        if r is not None and not J.contains_zero(r.Psi):
            nexc += 1
    out["rows"].append({"box_side": w, "boxes": len(bx),
                        "excluded_by_Psi": nexc,
                        "excluded_fraction": nexc / len(bx)})
    print("w=%.3f  boxes=%d  excluded=%.1f%%  (%.0fs)"
          % (w, len(bx), 100.0 * nexc / len(bx), time.time() - t0), flush=True)
    finest = (w, bx)

# ---- overestimation ratios and float sign map at the finest w --------------
w, bx = finest
bf = float(BETA)
ratios_v, ratios_g = [], []
sign_cells = 0
NS = 9                      # float samples per axis inside a box
for idx, (s0, s1, d0, d1) in enumerate(bx):
    vals = []
    grads = []
    pos = neg = False
    for i in range(NS):
        for k in range(NS):
            s_ = float(s0) + (float(s1) - float(s0)) * (i + 0.5) / NS
            d_ = float(d0) + (float(d1) - float(d0)) * (k + 0.5) / NS
            v = J.f_Dsep(bf, s_ + d_, s_)
            vals.append(v)
            if v > 0:
                pos = True
            elif v < 0:
                neg = True
            h = 1e-6
            grads.append((J.f_Dsep(bf, s_ + d_ + h, s_) - v) / h)
    if pos and neg:
        sign_cells += 1
    if idx % 7 == 0:        # subsample for the interval ratios
        T1 = iv.mpf([s0, s1])
        Dm = iv.mpf([d0, d1])
        T0 = T1 + Dm
        br = S.branches_for(B, T0, T1)
        r = S.evaluate(B, T0, T1, br)
        if r is None:
            continue
        tv = max(vals) - min(vals)
        tg = max(grads) - min(grads)
        if tv > 0:
            ratios_v.append(J.width(r.Dv) / tv)
        if tg > 0:
            ratios_g.append(J.width(r.D0) / tg)

ratios_v.sort()
ratios_g.sort()
out["finest"] = {
    "box_side": w,
    "float_samples_per_box": NS * NS,
    "median_enclosure_over_true_range_Dsep": ratios_v[len(ratios_v) // 2],
    "median_enclosure_over_true_range_dDsep_dth0": ratios_g[len(ratios_g) // 2],
    "ratio_sample_boxes": len(ratios_v),
    "cells_with_float_sign_change": sign_cells,
    "cells_total": len(bx),
    "fraction_cells_with_sign_change": sign_cells / len(bx),
}
out["seconds"] = time.time() - t0
with open(os.path.join(HERE, "conditioning.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out["finest"], indent=1))
print("total %.0fs" % out["seconds"])
