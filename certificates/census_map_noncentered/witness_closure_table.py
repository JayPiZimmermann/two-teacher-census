"""Per-witness CLOSURE table for the noncentered separated census.

The standing residual obligation of every face census (CERTIFICATE 9.2) is
that the separated half rests on a float locator: every family it returns is
certified, but a family it never seeds is not excluded, so a per-face census
was only ever "contains at least the rows listed".

A witness is CLOSED here when the sign-law-pruned branch-and-bound
(`face_bb_signlaw.py`) finishes with `complete = true` AND zero undecided
boxes: every point of the `(th1, D)` domain is then either excluded by a
certificate or inside a Krawczyk box, so the certified family list is
COMPLETE at that teacher -- exactly, not at least.  At such a witness the
obligation is discharged.

The table also checks the count against `faces_fixed.json`.  A mismatch is a
finding either way and is printed, never silently reconciled: the enumeration
and the locator are independent instruments.

WHAT THIS SUPPLIES TO THE COUNT-CONSTANCY SCHEMA.  `CountConstancyJ/
FaceCountSchemaJ.lean` proves face-wide count constancy from four inputs —
compact `K`, preconnected `S`, boundary containment (`hint`: every zero in
`K` lies in `interior K`), and fold-freeness (`hfold`) — plus a zero count
`N` at one witness (`hN`).  It explicitly does NOT assert the enumeration
interface, leaving `N` and its relation to the census family count to this
lane.  Two of those inputs are what a closed run here produces:

* `hN`.  The schema counts zeros of `censusAngleMapJ` in LABELLED student
  coordinates `x = (th0, th1)`, while this enumeration counts UNORDERED
  families on `(th1, D)`.  Each unordered pair `{th0, th1}` with
  `th0 != th1` gives two labelled zeros, so
      `N = 2 * (n_families_certified + n_exact_fit)`,
  reported below as `N_labelled`.  With exactly one exact-fit class at every
  witness so far this is `2*(n_separated + 1)`, the relation the schema
  records as measured — now backed by a CERTIFIED count rather than a float
  one.  The equivalence of the two zero sets is machine-checked:
  `SeparatedCount/CensusResidualJ.lean` identifies this enumeration's
  residual pair with `(generalJAngleEq0, generalJAngleEq1)`, which is
  `separatedAngleMapJ` by definition.
* `hint`.  A closed run has NO undecided box and no certified box touching
  the domain boundary — that is exactly what the seam offset, the `D`-range
  padding and the atom-line straddling were for, and it is why the run
  closes at all.  So boundary containment holds for THIS domain.

Two scope limits, both real.  The schema's `K` lives in `(th0, th1)` and
this domain in `(th1, D)`; that is a linear change of variables, but
transferring compactness and the interior statement across it is the schema
lane's step, not asserted here.  And `hfold` is not established here: the
per-family Schur data is recorded, and the tree's
`separatedAngleJacDetJ_ne_zero_of_kernelTeacherSchurDetT_ne_zero` is the
route from it, but this table does not run that argument.

Scope this table does NOT claim: closure at a witness is closure at a POINT.
Transporting it across a face needs the face to be count-homogeneous (F5 is
NOT -- see `f5_count_band.json`) plus either a witness per homogeneous
sub-region or a parametrized-teacher run over a `(beta, y)` box.

Usage: python3 witness_closure_table.py [prefix]   -> witness_closure.json
"""
import glob
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    pref = sys.argv[1] if len(sys.argv) > 1 else "L_"
    faces = json.load(open(os.path.join(HERE, "faces_fixed.json")))
    cens = {(round(w["beta"], 6), round(w["y"], 6)): w
            for w in faces["witnesses"]}
    rows = []
    for f in sorted(glob.glob(os.path.join(HERE,
                                           "signlaw_bb_%s*.json" % pref))):
        d = json.load(open(f))
        dd = d.get("dedup") or {}
        key = (round(d["beta"], 6), round(d["y"], 6))
        w = cens.get(key)
        closed = bool(d["complete"] and d["n_undecided"] == 0
                      and dd.get("disjointness_certified"))
        rows.append({
            "tag": os.path.basename(f),
            "beta": d["beta"], "y": d["y"],
            "face": w["face"] if w else "collar",
            "closed": closed,
            "complete": d["complete"],
            "n_undecided": d["n_undecided"],
            "undecided_area": d["undecided_area"],
            "steps": d["steps"],
            "seconds": d["seconds"],
            "n_families_certified": dd.get("n_separated_families"),
            "n_exact_fit": dd.get("n_exact_fit_classes"),
            "disjointness_certified": dd.get("disjointness_certified"),
            "n_families_locator": w["n_separated"] if w else None,
            "prefilter_area": d.get("prefilter_area"),
            "excluded_centered_form": d.get("excluded_centered_form"),
            "excluded_plain_form": d.get("excluded_plain_form"),
        })
        nf, ne = (rows[-1]["n_families_certified"],
                  rows[-1]["n_exact_fit"])
        rows[-1]["N_labelled"] = (2 * (nf + ne)
                                  if (nf is not None and ne is not None)
                                  else None)
        rows[-1]["count_agrees"] = (
            rows[-1]["n_families_certified"] == rows[-1]["n_families_locator"]
            if (w and rows[-1]["n_families_certified"] is not None) else None)

    print("%-6s %-7s %-8s %-7s %-6s %-9s %-7s %-10s %s"
          % ("face", "beta", "y", "closed", "fams", "locator", "agree",
             "N_labelled", "undecided"))
    for r in sorted(rows, key=lambda z: (z["face"], z["beta"], z["y"])):
        print("%-6s %-7.3f %-8.4f %-7s %-6s %-9s %-7s %-10s %d (%.1e)"
              % (r["face"], r["beta"], r["y"], r["closed"],
                 r["n_families_certified"], r["n_families_locator"],
                 r["count_agrees"], r["N_labelled"], r["n_undecided"],
                 r["undecided_area"]))

    per_face = {}
    for r in rows:
        per_face.setdefault(r["face"], []).append(r)
    summary = {f: {"n_witnesses": len(v),
                   "n_closed": sum(1 for z in v if z["closed"]),
                   "all_closed": all(z["closed"] for z in v),
                   "counts_certified": sorted(
                       {z["n_families_certified"] for z in v}),
                   "counts_locator": sorted({z["n_families_locator"]
                                             for z in v
                                             if z["n_families_locator"]
                                             is not None}),
                   "all_counts_agree": all(z["count_agrees"] for z in v
                                           if z["count_agrees"] is not None)}
              for f, v in per_face.items()}
    doc = {"object": "per-witness closure of the separated enumeration",
           "schema_interface": (
               "N_labelled = 2*(families + exact_fit) is the zero count of "
               "censusAngleMapJ in labelled coordinates, the hN input of "
               "CountConstancyJ/FaceCountSchemaJ; the zero-set identification "
               "is machine-checked in SeparatedCount/CensusResidualJ.lean"),
           "closed_means": "complete AND zero undecided boxes AND "
                           "disjointness of the certified classes certified",
           "n_witnesses": len(rows),
           "n_closed": sum(1 for r in rows if r["closed"]),
           "per_face": summary, "rows": rows}
    json.dump(doc, open(os.path.join(HERE, "witness_closure.json"), "w"),
              indent=1)
    print()
    print(json.dumps(summary, indent=1))
    print("CLOSED %d of %d witnesses" % (doc["n_closed"], len(rows)))
    mism = [r for r in rows if r["count_agrees"] is False]
    if mism:
        print("COUNT MISMATCHES (enumeration vs locator):")
        for r in mism:
            print("   beta=%.3f y=%.4f  certified=%s  locator=%s"
                  % (r["beta"], r["y"], r["n_families_certified"],
                     r["n_families_locator"]))


if __name__ == "__main__":
    main()
