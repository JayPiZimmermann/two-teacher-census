"""Merge every certified artefact of this lane into arrangement.json."""
import hashlib
import json
import os

def load(p, default=None):
    if not os.path.exists(p):
        return default
    with open(p) as fh:
        return json.load(fh)


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


faces = load("arrangement_faces.json")
collars = load("collars.json")
mosaic = load("mosaic_comparison.json")
band = load("band2.json")
quarter = {"locus": load("quarter_locus.json"),
           "saddle_certificate": load("quarter_saddle.json")}
sep_mixed = load("sep_mixed.json")
sep_fine = load("sep_fine.json")
sep_pos = load("sep_positive.json")
widget = load("face_widget.json")

if not faces or faces.get("face_count") != 10 or faces.get("distinct_census_values") != 4:
    raise RuntimeError("arrangement_faces.json must certify 10 faces and 4 signatures")
face_pairs = faces.get("visual_face_pairs") or {}
if (face_pairs.get("count") != 5
        or face_pairs.get("connectivity_checked") is not False
        or len(face_pairs.get("pairs") or []) != 5):
    raise RuntimeError("arrangement_faces.json must carry five non-topological visual pairs")
widget_pairs = ((widget or {}).get("orthogonal_column") or {}).get("visual_face_pairs") or {}
if (widget_pairs.get("count") != 5
        or widget_pairs.get("connectivity_checked") is not False
        or not widget_pairs.get("exact_pair_labels_agree")
        or len(widget_pairs.get("pairs") or []) != 5
        or not all(pair.get("check_passed") for pair in widget_pairs.get("pairs") or [])):
    raise RuntimeError("face_widget.json must pass all five local visual-pair witnesses")

input_files = [
    "arrangement_faces.json", "collars.json", "mosaic_comparison.json",
    "band2.json", "quarter_locus.json", "quarter_saddle.json",
    "sep_mixed.json", "sep_fine.json", "sep_positive.json",
    "face_widget.json",
]

out = {
    "title": "Certified planar arrangement and source/renderer census agreement "
             "for the CENTERED two-student map",
    "produced_by": "certificates/census_map_centered",
    "inputs": {
        "certificate": "certificates/census_map_centered/cert "
                       "(certificate.json, CERTIFICATE.md, centered.py)",
        "mosaic": "website/census-map.js (depth 8) and "
                  "website/precompute_census_map.js",
        "classifier": "website/widgets.js, loaded through the same DOM shim the "
                      "mosaic generator uses",
        "sha256": {
            "merged_artifacts": {p: sha256(p) for p in input_files
                                 if os.path.exists(p)},
            "source_files": (widget or {}).get("sources_sha256", {}),
            "producer": {"assemble.py": sha256("assemble.py")},
        }},
    "precision": "mpmath.iv, 200 bits with directed rounding (mp.prec 220) for "
                 "every certified statement; mpmath at 60 decimal digits for the "
                 "second-variation tests; float64 only as a locator, never as a "
                 "verdict",
    "validity_domain": "beta in [0, pi], y in R/2Z, centered kernel only "
                       "(phiCos, pi-periodic); two teacher atoms, two students",

    "headline": {
        "faces_of_the_arrangement": 10,
        "distinct_census_values_on_faces": 4,
        "explicit_non_topological_visual_face_pairs": 5,
        "finite_grid_flood_fill_groups_reported_by_the_shipped_mosaic":
            ((widget or {}).get("mosaic") or {}).get("flood_fill_groups"),
        "visual_pair_claim_scope":
            "presentation labels plus equal face signatures and the listed local "
            "exact beta=pi/2 witnesses; no connectivity, adjacency, or "
            "path-exclusion claim",
        "shipped_mosaic_state": "website/census-map.js at the SHA-256 recorded "
            "in inputs.sha256.source_files: 4 census keys in 5 finite-grid "
            "flood-fill groups; every renderer representative reproduces the "
            "certified face census where it lands (face_widget.json).  This is "
            "renderer agreement, not a topology certificate.",
        "history":
            "The asset of commit 72fb166 reported 8 finite-grid flood-fill "
            "groups: 4 cutoff artifacts and a coalesced visual grouping across "
            "the unresolved potential-curve pinch.  Both renderer defects were "
            "measured here "
            "(mosaic_comparison, sel_eps_artifact_band) and fixed in website/: "
            "the selector cutoff became scale-relative (commit c975b0e) and "
            "the generator takes the certified potential curve as a separator "
            "(commit f430fab).  The exact beta = pi/2 source census was later "
            "corrected to insert both common torque roots; it supplies local "
            "exact-column witnesses for five visual face pairs.  No "
            "connected-component conclusion is drawn from those witnesses."},

    "topology": (faces or {}).get("topology"),
    "visual_face_pairs": (faces or {}).get("visual_face_pairs"),
    "orthogonal_column_source_regression": {
        "artifact": "face_widget.json",
        "embedded_at": "widget_census_at_the_face_representatives.orthogonal_column",
        "cases_disagreeing": ((widget or {}).get("orthogonal_column") or {}).get(
            "cases_disagreeing"),
        "visual_face_pairs": ((widget or {}).get("orthogonal_column") or {}).get(
            "visual_face_pairs"),
    },
    "projective_beta_seam_source_regression": {
        "artifact": "face_widget.json",
        "embedded_at": "widget_census_at_the_face_representatives.beta_seam",
    },

    "why_the_census_is_a_function_of_the_three_determinants": {
        "statement":
            "Away from a common-root degenerate column, at a torque root t the "
            "mass direction is (s0,s1) = kappa * "
            "(-h(t-beta), h(t)), so P(t) = -kappa Wpot, W(t) = -kappa Wwgt, "
            "tau(t) = +kappa Wtau (the last using Wtau + Wpot = 2 Wwgt).  "
            "widgets.js calls a coincidence row a spurious minimum iff tau*P > 0 "
            "and W a b tau < 0, so kappa cancels and",
        "trap_present": "Wtau * Wpot < 0 at the root",
        "label_positive": "Wtau * Wwgt > 0",
        "label_mixed": "Wtau * Wwgt < 0",
        "consequence":
            "the census signature is a function of (number of torque roots; sign "
            "of each of Wtau, Wpot, Wwgt at each root).  The root count changes "
            "only where Wtau = 0 (the torque Wronskian is the Jacobian of "
            "t -> y, so its zero set is the fold locus of the covering), and the "
            "signs change only where one of the three vanishes.  That is exactly "
            "the certificate's boundary set off the degenerate columns.  At "
            "beta = pi/2 the mass-direction ratio is 0/0 and the two common "
            "roots are evaluated directly; see orthogonal_column_source_regression."},

    "faces": faces,
    "widget_census_at_the_face_representatives": widget,
    "collars": collars,
    "quarter_locus": quarter,
    "mosaic_comparison": mosaic and dict(mosaic, measured_against=
        "website/census-map.js as of commit 72fb166 (the 8-piece asset); "
        "both defects fixed in c975b0e and f430fab -- see headline.history"),
    "sel_eps_artifact_band": {
        "cutoff": "SEL_EPS = 1e-7, widgets.js coincidenceTypeAt, applied as an "
                  "ABSOLUTE cutoff to the factor W = s0|sin t| + s1|sin(t-beta)| "
                  "at the coincidence root",
        "measured_band_halfwidths_in_y": band,
        "mosaic_cell_height": 2.0 / 3072,
        "reading":
            "W is proportional to Wwgt, which vanishes exactly on y = 0 and "
            "y = +-1, so |W| < 1e-7 is a band around those two strata.  Its width "
            "is above one depth-8 cell only for beta < ~0.0125 and beta > ~3.129, "
            "which is exactly where the mosaic grew its four 514-cell slivers "
            "(measured extent beta in [0, 0.01636] and [3.12523, pi]).  The band "
            "is NONZERO at every beta (1e-6 to 1e-5 in the middle of the range), "
            "so the artifact is present everywhere and is merely unresolved: the "
            "agreement of the finite-grid group count at depths 6, 7 and 8 is "
            "not evidence of "
            "convergence, and a depth around 14-15 would resolve the band along "
            "the whole of y = 0 and y = +-1 and change the count."},
    "separated_families_mixed_sign": sep_mixed and {
        "stats": sep_mixed["stats"], "hits": sep_mixed["hits"],
        "grid": {"nb": sep_mixed["nb"], "ny": sep_mixed["ny"],
                 "dps": sep_mixed["dps"]},
        "verdict": "0 separated critical families with both student weights "
                   "nonzero were found for any of the 6144 mixed-sign teachers, "
                   "hence 0 local minima; since closed in Lean: "
                   "centered_mixedSign_distinct_critical_isExactFit proves every "
                   "such configuration is the exact fit, so this sweep is "
                   "corroboration, not the evidence",
        "fine_and_sensitivity_checks": sep_fine},
    "separated_families_positive_control": sep_pos and {
        "stats": sep_pos["stats"], "hits": sep_pos["hits"],
        "grid": {"nb": sep_pos["nb"], "ny": sep_pos["ny"], "dps": sep_pos["dps"]}},

    "containment_status": [
        "COINCIDENT families: the mass-free boundary equations are proved in "
        "Planar/GeneralTeachers/PhaseBoundary.lean and PotentialBoundary.lean "
        "(necessity and sufficiency, with the realising teacher).",
        "SAME-SIGN teachers (y < 0, both masses positive after the fold): "
        "centered_twoTeacher_generalMass_critical_iff_labeled_families "
        "(Planar/GeneralTeachers/Classification/FamilyLabels.lean) proves, for "
        "0 < beta < pi and positive teacher masses, with the STUDENT weights "
        "free: critical iff coincidence-line branch or quarter branch or "
        "separated branch, and the only family that can be a nonglobal local "
        "minimum is the coincidence line at a local maximum of the teacher "
        "potential with opposite-signed student weights "
        "(GenLineLocalMaxOuterSplit).",
        "MIXED-SIGN teachers (y > 0): "
        "centered_mixedSign_distinct_critical_isExactFit "
        "(Planar/GeneralTeachers/SeparatedRigidity/MixedSeparated.lean) proves "
        "that every critical configuration with distinct students mod pi and "
        "both weights nonzero is the exact fit; "
        "centered_mixedSign_no_quarter_critical_point rules out quarter-gap "
        "configurations (beta != pi/2).",
        "The quarter family on y = -1/2 (equal teacher masses): "
        "genQuarterBranch_eqTeacherMass_not_isLocalMin "
        "(Planar/GeneralTeachers/Benign/GenQuarter/EqTeacherMass.lean) proves "
        "it is never a local minimum, on ALL of (0, pi) \\ {pi/2} -- the "
        "22504-box interval certificate quarter_saddle.json is redundant "
        "corroboration on [0.002, pi - 0.002].",
        "The three critical beta values beta1*, pi/2, beta2* are not covered by "
        "any interval box (they cannot be: the deciding inequality's margin is "
        "O(distance)).  They are handled exactly in section 5 of the "
        "certificate.  The beta = pi/2 column's local-minimum census is "
        "additionally checked from its exact root formulas against the source "
        "classifier in face_widget.json.  The checker compares five explicit "
        "visual face-pair labels with the listed local exact-column witnesses and makes no "
        "connectivity or adjacency claim.",
    ],
}
with open("arrangement.json", "w") as fh:
    json.dump(out, fh, indent=1)
print("wrote arrangement.json (%d bytes)" % os.path.getsize("arrangement.json"))
