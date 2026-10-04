"""Certified-vs-shipped census scan across the ANTIPODAL BAND.

Usage: python3 band_scan.py <slot> <nslots> <b0> <b1> <nb> <ny> [tag]
Writes band<tag>_<slot>.json
"""
import json
import os
import sys
import time

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
slot, nslots = int(sys.argv[1]), int(sys.argv[2])
B0, B1 = float(sys.argv[3]), float(sys.argv[4])
NB, NY = int(sys.argv[5]), int(sys.argv[6])
TAG = sys.argv[7] if len(sys.argv) > 7 else ""

pts = []
for i in range(NB):
    beta = B0 + (B1 - B0) * (i + 0.5) / NB
    for k in range(NY):
        y = -1 + 2 * (k + 0.5) / NY
        pts.append((round(beta, 10), round(y, 10)))
mine = [p for n, p in enumerate(pts) if n % nslots == slot]

out = []
t0 = time.time()
for beta, y in mine:
    try:
        r = X.certified_census(beta, y, N=180)
    except Exception as e:                      # noqa
        out.append({"beta": beta, "y": y, "error": str(e)})
    else:
        out.append({
            "beta": beta, "y": y, "census": r["census"],
            "n_torque_roots": r["n_torque_roots"],
            "torque_t": [q["t_mid"] for q in r["torque_roots"]],
            "n_separated": r["n_separated"],
            "n_spurious_sep": sum(1 for s in r["separated"]
                                  if s["schur"] == "spurious"),
            "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
            "sep_schur": [s["schur"] for s in r["separated"]],
            "certified": r["certified_modulo_locator_completeness"],
            "torque_undecided": len(r["torque_undecided"]),
            "locator_failed": len(r["locator_failed_certification"]),
        })
    with open(os.path.join(HERE, "band%s_%d.json" % (TAG, slot)), "w") as fh:
        json.dump(out, fh)
print("slot %d: %d points, %.0fs" % (slot, len(out), time.time() - t0))
