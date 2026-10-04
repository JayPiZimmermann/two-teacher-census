"""Certified spot check of the two outer beta intervals (0, beta*_1) and
(beta*_2, 2pi), whose adaptive coverings did not complete.

For each centre, certify_beta_box is run on a shrinking box about it, and the
WIDEST box that certifies is recorded together with its certified root-count
pair (Wtau_total, Wpot_total).  This is a SPOT CHECK, not a covering: it
certifies the counts on the accepted boxes only, and their total measure is
reported against the combined length of the two intervals.

Writes spot_check.json.
"""
import json
import os
import time

import mpmath
from mpmath import iv

import noncentered as J
import certify_coincident as C

HERE = os.path.dirname(os.path.abspath(__file__))
J.set_prec(C.PREC)

u, B1, B2 = C.certify_beta_star()
b1 = float(J.mid(B1))
b2 = float(J.mid(B2))
two_pi = float(J.mid(2 * J.PI_IV()))

CENTRES = [0.01, 0.05, 0.10, 0.30, 0.60, 1.00, 1.50, 2.00, 2.20,
           4.10, 4.60, 5.00, 5.60, 6.10, 6.24, 6.27]

out = {"intervals": {"(0, beta*_1)": [0.0, b1], "(beta*_2, 2pi)": [b2, two_pi]},
       "interval_total_length": (b1 - 0.0) + (two_pi - b2),
       "rows": []}
t00 = time.time()
for c in CENTRES:
    crits = [0.0, b1, b2, two_pi, float(J.mid(J.PI_IV()))]
    dist = min(abs(c - z) for z in crits)
    # cap the first box: certify_beta_box's sinc series is guarded at |z| <= 3.4,
    # and a very wide beta box pushes the reduced-form arguments past it
    hw = min(dist / 2, 0.25)
    accepted = None
    counts = None
    tries = 0
    while hw > dist / 2 ** 19:
        tries += 1
        B = iv.mpf([mpmath.mpf(repr(c)) - mpmath.mpf(repr(hw)),
                    mpmath.mpf(repr(c)) + mpmath.mpf(repr(hw))])
        r = C.certify_beta_box(B)
        if r.get("ok"):
            accepted = hw
            counts = (r["Wtau_total"], r["Wpot_total"])
            break
        hw = hw / 2
    out["rows"].append({"centre": c, "certified_half_width": accepted,
                        "counts_Wtau_Wpot": counts, "tries": tries})
    print("centre %-5s  hw=%s  counts=%s" % (c, accepted, counts), flush=True)

ok_rows = [r for r in out["rows"] if r["certified_half_width"]]
out["n_certified"] = len(ok_rows)
out["n_centres"] = len(CENTRES)
out["certified_measure"] = sum(2 * r["certified_half_width"] for r in ok_rows)
out["fraction_of_intervals"] = out["certified_measure"] / out["interval_total_length"]
out["seconds"] = time.time() - t00
with open(os.path.join(HERE, "spot_check.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print("certified %d/%d centres, measure %.4f of %.4f (%.2f%%), %.0fs"
      % (out["n_certified"], out["n_centres"], out["certified_measure"],
         out["interval_total_length"], 100 * out["fraction_of_intervals"],
         out["seconds"]))
