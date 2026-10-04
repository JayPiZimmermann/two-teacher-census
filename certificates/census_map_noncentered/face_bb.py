"""Per-teacher CERTIFIED-COMPLETE census at the seven face representatives and
their witnesses: the separated stratum enumerated by the two-equation interval
branch and bound of census_cert.separated_families (exclusion on the
radial-solved residual pair (G0, G1), Krawczyk existence/uniqueness on the
survivors) -- NOT by the float locator.  Closing this discharges the locator-
completeness obligation at every point it certifies: the census then lists
EXACTLY the separated families that exist outside the coincidence collar
D < delta.

Usage: python3 face_bb.py <slot> <nslots>     -> face_bb_<slot>.json
"""
import json
import os
import sys
import time

import mpmath

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
FACES = [
    ("F1", (0.75, 0.40)), ("F1", (1.50, 0.60)), ("F1", (2.60, 0.20)),
    ("F1", (0.50, 0.80)),
    ("F2", (0.75, -0.02)), ("F2", (0.30, -0.02)), ("F2", (1.50, -0.04)),
    ("F2", (2.20, -0.06)),
    ("F3", (3.10, -0.06)), ("F3", (3.10, -0.10)), ("F3", (3.13, -0.08)),
    ("F3", (3.05, -0.09)),
    ("F4", (0.75, -0.40)), ("F4", (0.30, -0.50)), ("F4", (1.50, -0.60)),
    ("F4", (2.00, -0.30)),
    ("F5", (2.60, -0.50)), ("F5", (2.80, -0.50)), ("F5", (3.00, -0.40)),
    ("F5", (2.40, -0.50)),
    ("F6", (0.75, -0.98)), ("F6", (0.30, -0.98)), ("F6", (1.50, -0.96)),
    ("F6", (2.20, -0.94)),
    ("F7", (3.10, -0.94)), ("F7", (3.10, -0.90)), ("F7", (3.13, -0.92)),
    ("F7", (3.05, -0.91)),
]

slot, nslots = int(sys.argv[1]), int(sys.argv[2])
BUDGET = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
DELTA = mpmath.mpf("0.02")
mine = [j for n, j in enumerate(FACES) if n % nslots == slot]
out = []
for fid, (b, y) in mine:
    t0 = time.time()
    try:
        r = X.census_of(b, y, sep=True, delta=DELTA, budget=BUDGET)
        out.append({
            "face": fid, "beta": b, "y": y,
            "census": r["census"], "certified": r["certified"],
            "n_torque_roots": r["n_torque_roots"],
            "n_separated": r["n_separated"],
            "families": [{"th1": s["th1_mid"], "D": s["D_mid"],
                          "schur": s["schur"], "exact_fit": s.get("exact_fit"),
                          "sgn_c": [s["sgn_c0"], s["sgn_c1"]],
                          "width": s["enclosure_width"]}
                         for s in r["separated"]],
            "n_undecided": len(r["separated_undecided"]),
            "undecided": r["separated_undecided"][:40],
            "undecided_area": sum(u["area"] for u in r["separated_undecided"]),
            "steps": r["separated_steps"],
            "collar_delta": float(DELTA),
            "collar_area": r["separated_collar_area"],
            "seconds": time.time() - t0,
        })
    except Exception as e:                     # noqa
        out.append({"face": fid, "beta": b, "y": y, "error": repr(e),
                    "seconds": time.time() - t0})
    with open(os.path.join(HERE, "face_bb_%d.json" % slot), "w") as fh:
        json.dump(out, fh, indent=1)
    print("%s (%.2f, %+.2f): %s  undec=%d  %.0fs"
          % (fid, b, y, out[-1].get("census", "ERROR"),
             out[-1].get("n_undecided", -1), out[-1]["seconds"]), flush=True)
print("slot %d done %d" % (slot, len(out)))
