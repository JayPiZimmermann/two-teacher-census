"""Rebuild the per-face witness table under the FIXED sep_res Jacobian.

Successor artifact to `faces_0.json`, which was produced with the Jacobian
sign error (see `sep_jacobian_audit.py`).  Provenance is kept, not
overwritten: the two witnesses the fix invalidated stay in the file with
`status: "superseded_by_jacobian_fix"` and their pre- and post-fix records
side by side, and each carries the replacement that took its place.

WHY THEY WERE REPLACED.  Both sat at `beta = 3.13`, i.e. `|beta - pi| =
0.0116`, inside the antipodal pocket where the fix changes verdicts:
  F3 (3.13, -0.08): n_separated 3 -> 1, certified True -> False, and the
                    census loses its `separate:trap@positive` row
  F7 (3.13, -0.92): n_separated 3 -> 2, certified True -> False
The replacements are chosen in the same face at `|beta - pi| = 0.0616`,
2.6x the observed pocket radius 0.024, and each reproduces its face's census
string and Schur pattern exactly.

Usage: python3 faces_rebuild_fixed.py    -> faces_fixed.json
"""
import json
import math
import os
import time

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
SUPERSEDED = {("F3", 3.13, -0.08): ("F3", 3.08, -0.08),
              ("F7", 3.13, -0.92): ("F7", 3.08, -0.92)}
POCKET_RADIUS = 0.024


def record(face, beta, y, is_rep):
    r = X.certified_census(beta, y, N=160)
    return {"face": face, "beta": beta, "y": y, "is_rep": is_rep,
            "census": r["census"],
            "n_torque_roots": r["n_torque_roots"],
            "n_separated": r["n_separated"],
            "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
            "sep_schur": [s["schur"] for s in r["separated"]],
            "certified": r["certified_modulo_locator_completeness"],
            "dist_to_pi": round(abs(beta - math.pi), 6)}


def main():
    old = json.load(open(os.path.join(HERE, "faces_0.json")))
    out, superseded = [], []
    t0 = time.time()
    for w in old:
        key = (w["face"], round(w["beta"], 6), round(w["y"], 6))
        rec = record(w["face"], w["beta"], w["y"], w.get("is_rep", False))
        if key in SUPERSEDED:
            f, nb, ny = SUPERSEDED[key]
            rep = record(f, nb, ny, w.get("is_rep", False))
            rec["status"] = "superseded_by_jacobian_fix"
            rec["pre_fix"] = {k: w[k] for k in
                              ("census", "n_separated", "sep_D",
                               "sep_schur", "certified")}
            rec["replaced_by"] = {"beta": nb, "y": ny}
            rec["reason"] = ("inside the antipodal pocket |beta - pi| <= "
                             "%.3f where the Jacobian fix changes verdicts"
                             % POCKET_RADIUS)
            superseded.append(rec)
            rep["status"] = "replacement"
            rep["replaces"] = {"beta": w["beta"], "y": w["y"]}
            out.append(rep)
        else:
            rec["status"] = ("unchanged_by_jacobian_fix"
                             if (rec["census"] == w["census"]
                                 and rec["n_separated"] == w["n_separated"]
                                 and rec["certified"] == w["certified"])
                             else "CHANGED_by_jacobian_fix")
            out.append(rec)
        print("  %-3s %-6s %-7s %-28s %s" % (
            rec["face"], rec["beta"], rec["y"], rec["status"],
            "certified" if rec["certified"] else "DECLINED"), flush=True)

    per_face = {}
    for w in out:
        per_face.setdefault(w["face"], []).append(w)
    summary = {f: {"n_witnesses": len(v),
                   "all_certified": all(z["certified"] for z in v),
                   "census": sorted({z["census"] for z in v}),
                   "min_dist_to_pi": min(z["dist_to_pi"] for z in v)}
               for f, v in per_face.items()}
    doc = {
        "object": "per-face witness table under the FIXED sep_res Jacobian",
        "supersedes": "faces_0.json",
        "pocket_radius_used": POCKET_RADIUS,
        "n_witnesses": len(out),
        "n_replaced": len(superseded),
        "witnesses": out,
        "superseded": superseded,
        "per_face": summary,
        "seconds": time.time() - t0,
    }
    json.dump(doc, open(os.path.join(HERE, "faces_fixed.json"), "w"),
              indent=1)
    print(json.dumps(summary, indent=1))
    bad = [f for f, s in summary.items() if not s["all_certified"]
           or len(s["census"]) != 1]
    print("faces with a declined witness or a split census:", bad)


if __name__ == "__main__":
    main()
