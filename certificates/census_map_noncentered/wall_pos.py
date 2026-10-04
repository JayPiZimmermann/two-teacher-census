"""Certify selected pointwise census witnesses in the mixed collar.

The locator-fixed classifier (widgets.js, commit 2201dd0) shows, at small beta
near the strata y = 0 and |y| = 1, a census with ONE separate:trap@mixed row
where sampled witnesses carrying the proposed F1 label have two. This scan
certifies that census value at
explicit points (locate-then-certify, geometric D-grid from 1e-4, every family
Krawczyk-certified) and records, per beta, adjacent sampled points with the
one-trap and two-trap census values. The samples do not certify an intervening
wall, its uniqueness, or its continuation in beta.

The historical filename ``wall_pos`` is retained for provenance; neither the
producer nor its JSON certifies a wall. The authoritative scope ledger is
``CERTIFICATE.md`` (status 2026-08-24).

Writes wall_pos.json.
"""
import json
import os
import sys
import time

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
BETAS = [0.15, 0.30, 0.50, 0.75, 1.00]
YS = [0.004, 0.008, 0.012, 0.020, 0.030, 0.050, 0.080, 0.120, 0.200]

slot, nslots = int(sys.argv[1]), int(sys.argv[2])
pts = [(b, y) for b in BETAS for y in YS]
mine = [p for n, p in enumerate(pts) if n % nslots == slot]
out = []
for beta, y in mine:
    t0 = time.time()
    try:
        r = X.certified_census(beta, y, N=260, dmin=1e-4)
        out.append({
            "beta": beta, "y": y, "census": r["census"],
            "n_torque_roots": r["n_torque_roots"],
            "n_separated": r["n_separated"],
            "n_spurious_sep": sum(1 for s in r["separated"]
                                  if s["schur"] == "spurious"
                                  and not s.get("exact_fit")),
            "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
            "sep_schur": [s["schur"] for s in r["separated"]],
            "certified": r["certified_modulo_locator_completeness"],
            "locator_failed": len(r["locator_failed_certification"]),
            "seconds": time.time() - t0,
        })
    except Exception as e:                      # noqa
        out.append({"beta": beta, "y": y, "error": repr(e)})
    with open(os.path.join(HERE, "wall_pos_%d.json" % slot), "w") as fh:
        json.dump(out, fh, indent=1)
    q = out[-1]
    print("(%.2f, %.3f): %s  [cert=%s, %.0fs]"
          % (beta, y, q.get("census", "ERR"), q.get("certified"),
             q.get("seconds", 0)), flush=True)
print("slot %d done" % slot)
