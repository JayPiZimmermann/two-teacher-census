"""Certified y-extent of the separated wall near y = 0 and y = +-1."""
import json, os, sys, time
import census_cert as X
HERE = os.path.dirname(os.path.abspath(__file__))
slot, nslots = int(sys.argv[1]), int(sys.argv[2])
BETAS = [0.30, 0.75, 1.50, 2.20, 2.60, 3.00, 3.10, 3.13]
YS = [-0.5, -0.3, -0.2, -0.15, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02,
      -0.01, -0.005, -0.002]
pts = [(b, y) for b in BETAS for y in YS]
mine = [p for n, p in enumerate(pts) if n % nslots == slot]
out = []
for beta, y in mine:
    try:
        r = X.certified_census(beta, y, N=160)
        out.append({"beta": beta, "y": y, "census": r["census"],
                    "n_sep": r["n_separated"],
                    "n_spur": sum(1 for s in r["separated"] if s["schur"] == "spurious"),
                    "sepD": [round(s["D_mid"], 6) for s in r["separated"]],
                    "certified": r["certified_modulo_locator_completeness"]})
    except Exception as e:
        out.append({"beta": beta, "y": y, "error": str(e)})
    with open(os.path.join(HERE, "ywall_%d.json" % slot), "w") as fh:
        json.dump(out, fh)
print("slot %d done %d" % (slot, len(out)))
