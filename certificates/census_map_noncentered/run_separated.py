"""Driver for TASK 3.

Usage: python3 run_separated.py <slot> <nslots> <halfwidth> <delta> <minw> <budget>
Certifies a list of narrow beta WINDOWS (not a covering of the beta axis --
that is stated explicitly in the certificate) and writes separated_<slot>.json.
"""
import json
import os
import sys
import time

import mpmath
from mpmath import iv

import noncentered as J
import certify_separated as S

HERE = os.path.dirname(os.path.abspath(__file__))

# window centres, spread over (0, 2pi), avoiding 0, pi, 2pi
CENTRES = ["0.30", "0.70", "1.10", "1.55", "1.95", "2.40", "2.80",
           "3.50", "3.90", "4.35", "4.80", "5.25", "5.70", "6.05"]


def main():
    slot = int(sys.argv[1])
    nslots = int(sys.argv[2])
    hw = mpmath.mpf(sys.argv[3])
    delta = mpmath.mpf(sys.argv[4])
    minw = mpmath.mpf(sys.argv[5])
    budget = int(sys.argv[6])
    out = []
    for i, c in enumerate(CENTRES):
        if i % nslots != slot:
            continue
        cm = mpmath.mpf(c)
        B = iv.mpf([cm - hw, cm + hw])
        t0 = time.time()
        r = S.halfdomain_bb(B, delta, minw, budget=budget)
        r["seconds"] = time.time() - t0
        r["centre"] = c
        out.append(r)
        print("beta %s +-%s : closed %.4f of the half domain, %d branch boxes, "
              "%d unresolved, %.0fs"
              % (c, mpmath.nstr(hw, 5), r["closed_fraction_of_half_domain"],
                 r["n_branch_boxes"], r["n_unresolved"], r["seconds"]),
              flush=True)
        with open(os.path.join(HERE, "separated_%d.json" % slot), "w") as fh:
            json.dump(out, fh)


if __name__ == "__main__":
    main()
