"""Certified y-extent of the separated wall in the SAME-SIGN sector y > 0."""
import json
import os
import sys

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
slot, nslots = int(sys.argv[1]), int(sys.argv[2])
BETAS = [0.15, 0.50, 1.00, 1.50, 2.20, 2.60, 3.00]
YS = [0.002, 0.005, 0.010, 0.020, 0.035, 0.050, 0.080, 0.120, 0.200, 0.350]
pts = [(b, y) for b in BETAS for y in YS]
mine = [p for n, p in enumerate(pts) if n % nslots == slot]
out = []
for beta, y in mine:
    try:
        r = X.certified_census(beta, y, N=160)
        out.append({"beta": beta, "y": y, "census": r["census"],
                    "n_sep": r["n_separated"],
                    "n_spur": sum(1 for s in r["separated"]
                                  if s["schur"] == "spurious"),
                    "sepD": [round(s["D_mid"], 6) for s in r["separated"]],
                    "certified": r["certified_modulo_locator_completeness"]})
    except Exception as e:                      # noqa
        out.append({"beta": beta, "y": y, "error": str(e)})
    with open(os.path.join(HERE, "ywallp_%d.json" % slot), "w") as fh:
        json.dump(out, fh)
print("slot %d done %d" % (slot, len(out)))
