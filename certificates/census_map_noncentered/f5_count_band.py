"""Sample separated-family counts near the proposed F5 label.

The finite grid contains named points with counts 3/3/3/1
(`faces_fixed.json`) and adjacent sampled pairs whose counts differ. Those
pairs are candidate-transition brackets only: they do not certify intervening
walls, a count-homogeneous global region, a required subdivision, or
constancy on any other proposed label.

At the sampled points carrying the proposed F5 label, the extra returned
families are saddles, so the
recorded trap-only census strings agree there despite different family counts.
This pointwise observation does not prove F5-wide census constancy or a global
census arrangement.

The rectangular grid crosses the sampled torque-lens classification. Points
with one `coincident:trap@mixed` row carry the proposed F4 label and are not
included in the candidate-transition brackets for the proposed F5 label. At
`beta=2.40`, this excludes the sampled rows at `y<=-0.54`.

Usage: python3 f5_count_band.py   -> f5_count_band.json

The historical filename is retained for provenance. ``CERTIFICATE.md``
(status 2026-08-24) is the authoritative scope ledger.
"""
import json
import os

import census_cert as X

HERE = os.path.dirname(os.path.abspath(__file__))
BETAS = [2.40, 2.60, 2.80, 3.00, 3.10]
YS = [-0.60, -0.58, -0.56, -0.54, -0.52, -0.50, -0.48, -0.46, -0.45,
      -0.44, -0.42, -0.40]
F5_COINCIDENT_TRAPS = 2


def main():
    rows = []
    for b in BETAS:
        for y in YS:
            r = X.certified_census(b, y, N=140)
            cen = r["census"]
            rows.append({
                "beta": b, "y": y, "census": cen,
                "n_coincident_traps": cen.count("coincident:trap"),
                "has_proposed_F5_label":
                    cen.count("coincident:trap") == F5_COINCIDENT_TRAPS,
                "n_separated": r["n_separated"],
                "schur": sorted({s["schur"] for s in r["separated"]}),
                "sep_D": [round(s["D_mid"], 6) for s in r["separated"]],
                "certified": r["certified_modulo_locator_completeness"]})
            print("  beta=%.2f y=%-7s proposed_F5=%-5s n_sep=%d schur=%s"
                  % (b, y, rows[-1]["has_proposed_F5_label"],
                     rows[-1]["n_separated"],
                     ",".join(rows[-1]["schur"])), flush=True)

    brackets = {}
    for b in BETAS:
        seq = [r for r in rows
               if r["beta"] == b and r["has_proposed_F5_label"]]
        seq.sort(key=lambda z: z["y"])
        lower = upper = None
        for i in range(len(seq) - 1):
            a, c = seq[i], seq[i + 1]
            if a["n_separated"] == 1 and c["n_separated"] == 3:
                lower = [a["y"], c["y"]]
            if a["n_separated"] == 3 and c["n_separated"] == 1:
                upper = [a["y"], c["y"]]
        band = [r["y"] for r in seq if r["n_separated"] == 3]
        brackets[str(b)] = {
            "lower_count_change_bracket": lower,
            "upper_count_change_bracket": upper,
            "count3_y_scanned": band,
            "sampled_count3_y_span": (round(max(band) - min(band), 4)
                                      if len(band) > 1 else 0.0),
            "n_sampled_points_with_proposed_F5_label": len(seq)}
        print("beta=%.2f: candidate lower %s  candidate upper %s  "
              "sampled count-3 y at %s"
              % (b, lower, upper, band), flush=True)

    schur_all = sorted({t for r in rows if r["has_proposed_F5_label"]
                        for t in r["schur"]})
    doc = {
        "object": "finite proposed-F5 family-count samples and adjacent "
                  "candidate-transition brackets; no global walls, inter-sample "
                  "continuation, or F5-wide constancy are certified",
        "scope_note": "Every row is pointwise. Proposed-label membership is a "
                      "sample classifier, not a certified global cell.",
        "authoritative_ledger": {
            "path": "CERTIFICATE.md",
            "status_date": "2026-08-24",
        },
        "proposed_F5_sample_identified_by":
            "two coincident:trap@mixed rows (one row = proposed F4, outside "
            "the torque lens)",
        "all_schur_types_seen_at_sampled_proposed_F5_points": schur_all,
        "census_string_constant_on_sampled_proposed_F5_points": len({
            r["census"] for r in rows
            if r["has_proposed_F5_label"]}) == 1,
        "candidate_transition_brackets": brackets,
        "rows": rows}
    json.dump(doc, open(os.path.join(HERE, "f5_count_band.json"), "w"),
              indent=1)
    print("\nall Schur types seen at sampled proposed-F5 points:", schur_all)
    print("census string constant on sampled proposed-F5 points:",
          doc["census_string_constant_on_sampled_proposed_F5_points"])


if __name__ == "__main__":
    main()
