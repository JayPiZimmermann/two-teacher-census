"""Certified census on a grid over the FUNDAMENTAL DOMAIN beta in (0, pi),
y in (-1, 1] of the noncentered map.  (The census is invariant under
beta -> 2pi - beta; see CERTIFICATE.md section on the symmetry.)

Usage: python3 map_scan.py <slot> <nslots>
"""
import json
import os
import sys
import time

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
slot, nslots = int(sys.argv[1]), int(sys.argv[2])

BETAS = [0.05, 0.15, 0.30, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00,
         2.20, 2.40, 2.60, 2.80, 3.00, 3.10]
YS = [-0.98, -0.95, -0.90, -0.80, -0.60, -0.40, -0.20, -0.12, -0.06, -0.02,
      0.02, 0.06, 0.12, 0.20, 0.40, 0.60, 0.80, 0.90, 0.95, 0.98, 1.00]

pts = [(b, y) for b in BETAS for y in YS]
mine = [p for n, p in enumerate(pts) if n % nslots == slot]
out = []
t0 = time.time()
for beta, y in mine:
    try:
        r = X.certified_census(beta, y, N=160)
        out.append({
            "beta": beta, "y": y, "census": r["census"],
            "n_torque_roots": r["n_torque_roots"],
            "n_separated": r["n_separated"],
            "n_spurious_sep": sum(1 for s in r["separated"]
                                  if s["schur"] == "spurious"),
            "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
            "sep_schur": [s["schur"] for s in r["separated"]],
            "sep_sgn_c": [[s["sgn_c0"], s["sgn_c1"]] for s in r["separated"]],
            "certified": r["certified_modulo_locator_completeness"],
        })
    except Exception as e:                     # noqa
        out.append({"beta": beta, "y": y, "error": str(e)})
    with open(os.path.join(HERE, "map_%d.json" % slot), "w") as fh:
        json.dump(out, fh)
print("slot %d: %d points, %.0fs" % (slot, len(out), time.time() - t0))
