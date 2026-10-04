"""
COUNT-MATCHING record for the separated-family completeness residual
(CERTIFICATE.md section 9.2).

The route: a Lean-proved upper bound B on the number of separated critical
families of a noncentered teacher closes locator completeness at every
witness where the locator FINDS B families (found >= bound => found = all).

This script is the wiring of that route.  It
  1. reads the proven per-sector bounds from PROVEN_BOUNDS below (each entry
     must name the Lean theorem it transcribes; None = not proved),
  2. re-runs the locator (numpy mirror of census_cert.locate_separated,
     validated in count_scan.py) at the 28 certified face witnesses of
     faces_0.json plus the collar witnesses,
  3. reads the per-face count maxima of the slack scan
     (separated_count_scan.json, skill-1t artifact),
  4. records per-face closure status in count_match.json.

SECTOR NAMING FIX (2026-08-08, third pass).  Earlier passes of this file
called `0 < y < 1` the "same_sign" sector and `-1 < y < 0` "mixed".  That is
BACKWARDS: `masses_at(y)` is `psi = (y+1)pi/2, (s0,s1) = (sin psi, cos psi)`,
so `y in (-1,0)` gives `psi in (0,pi/2)` and BOTH masses positive, while
`y in (0,1)` gives `psi in (pi/2,pi)` and `s0 > 0 > s1`.  The keys are now
`positive_masses` (y < 0: faces F2-F7) and `mixed_signs` (y > 0: face F1 and
the collar).  The per-face bound VALUES are unchanged; only the labels were
wrong.  This matters because the two sectors are governed by DIFFERENT Lean
theorems (straddle vs same-side, below).

Current state (2026-08-08, second pass): the strict SIGN LAW of the scalar
is PROVED in Lean (Planar/SignedN2/FreeMassJ/SeparatedCount/SignLawJ.lean:
`separatedNumJ_neg_of_offset_lt_gap`, `separatedNumJ_pos_of_gap_lt_offset`,
`separatedNumJ_eq_zero_iff_lattice`, all on the standard axiom set), with
first counting consequence `separatedPairCritical_student0/1_offsets_straddle`
(LatticeStraddleJ.lean): each student's two teacher-offsets straddle the
gap lattice at every positive-teacher separated family.

NO count bound is proved yet.  The remaining blocker for the per-sector
caps (tight values 4 same-sign / 3 mixed) is now NAMED PRECISELY
(measurements + kill-records in branch_cap_scan.py, 2026-08-08 second
pass): the BRANCH CAP is the all-line crossing bound of the load curve —
for every teacher, every D, and every direction (alpha, gamma) != 0, the
combination alpha*(s0*W(t,D)+s1*W(t-beta,D)) + gamma*(s0*h(t)+s1*h(t-beta))
has at most 4 sign changes per period (the student load is one line;
measured tight, and a global property of the curve).  Routes REFUTED by
slack computation (do not re-try): Sturm-vs-forcing (forcing has 8 sign
changes at 24+/28 witnesses => cap 10, slack 6), per-h_load-arc ratio
monotonicity (4 zeros bunch in one arc), per-arc disconjugacy (cap 24),
monotone direction / naive winding bookkeeping (load Wronskian has up to
6 sign changes; net winding differs by sector while the cap does not).
Beyond the branch cap the family caps also need the D-profile crossing
bound (fold-profile architecture).  Until those land, every face records
"open".
"""
import json
import math

import count_scan as CS

# Per-sector proven upper bounds on the number of separated critical
# families (GeneralJSeparatedPairCritical solutions modulo gauge, exact fit
# excluded).  Each value, when set, MUST be the transcription of a Lean
# theorem named here.  None = no such theorem exists in the tree.
PROVEN_BOUNDS = {
    # -1 < y < 0  =>  psi in (0, pi/2)  =>  s0 > 0 AND s1 > 0.
    # Faces F2-F7.  Needed tight value 3 (F4: 1).
    # Sign-law geometry available in Lean: each student's two teacher-offsets
    # STRADDLE the gap lattice (SeparatedCount/LatticeStraddleJ.lean,
    # separatedPairCritical_student0/1_offsets_straddle).
    "positive_masses": {"bound": None, "lean_theorem": None},
    #  0 < y < 1   =>  psi in (pi/2, pi) =>  s0 > 0 > s1.
    # Face F1 and the same-sign collar.  Needed tight value 4.
    # Sign-law geometry available in Lean: each student's two teacher-offsets
    # lie on the SAME open side (SeparatedCount/MixedSameSideJ.lean,
    # separatedPairCritical_student0/1_offsets_same_side).
    "mixed_signs": {"bound": None, "lean_theorem": None},
}


def sector_of(y):
    """Which mass sector the census parameter y lands in.  NOTE the sign
    convention: (s0,s1) = (sin psi, cos psi) with psi = (y+1)pi/2, so y < 0
    is the POSITIVE-mass sector and y > 0 the mixed-sign one."""
    return "mixed_signs" if 0 < y < 1 else "positive_masses"


def main():
    faces = json.load(open("faces_0.json"))
    try:
        scan = json.load(open("separated_count_scan.json"))
    except FileNotFoundError:
        scan = None

    witnesses = [(w["face"], w["beta"], w["y"], w["n_separated"])
                 for w in faces]
    witnesses += [("collar", 0.15, 0.004, 1), ("collar", 0.15, 0.996, 1)]

    rows = []
    for face, beta, y, n_cert in witnesses:
        fams = CS.count_families(beta, y, N=260)
        fams2 = CS.count_families(beta, y, N=520)
        sec = sector_of(y)
        pb = PROVEN_BOUNDS[sec]["bound"]
        closed = (pb is not None and len(fams2) == pb
                  and len(fams) == len(fams2))
        rows.append({
            "face": face, "beta": beta, "y": y,
            "sector": sec,
            "n_certified_faces0": n_cert,
            "n_locator": len(fams), "n_locator_2N": len(fams2),
            "proven_bound": pb,
            "closed": bool(closed),
        })
        print(face, beta, y, "found", len(fams), len(fams2),
              "bound", pb, "closed", closed, flush=True)

    per_face = {}
    for r in rows:
        per_face.setdefault(r["face"], []).append(r)
    summary = {}
    for face, rs in per_face.items():
        ns = [r["n_locator_2N"] for r in rs]
        scan_max = None
        if scan and scan.get("scan") and face in scan["scan"]:
            scan_max = scan["scan"][face]["max"]
        summary[face] = {
            "witness_counts": ns,
            "scan_max": scan_max,
            "needed_tight_bound": max([n for n in ns]
                                      + ([scan_max] if scan_max else [])),
            "status": ("closed" if all(r["closed"] for r in rs) else "open"),
        }

    out = {
        "proven_bounds": PROVEN_BOUNDS,
        "witnesses": rows,
        "per_face": summary,
        "sign_law": (
            "PROVED (2026-08-08): SeparatedCount/SignLawJ.lean, "
            "separatedNumJ_neg_of_offset_lt_gap / "
            "separatedNumJ_pos_of_gap_lt_offset / "
            "separatedNumJ_eq_zero_iff_lattice, axiom-clean"),
        "sign_law_geometry": (
            "BOTH sectors now covered in Lean: positive masses -> STRADDLE "
            "(LatticeStraddleJ.lean, 2026-08-08), mixed signs -> SAME SIDE "
            "(MixedSameSideJ.lean, 2026-08-08: "
            "separatedNumJ_products_nonneg_of_separatedPairCritical / "
            "separatedNumJ_no_opposite_open_sides / "
            "separatedPairCritical_student0/1_offsets_same_side / "
            "separatedNumJ_products_nonneg_nonvacuous, axiom-clean).  These "
            "license the ORDER-ONLY box exclusion measured in "
            "face_bb_signlaw.py.  They do NOT bound the family count: the "
            "straddle/same-side condition is IMPLIED by the load equation "
            "s0*num(D,a)+s1*num(D,b)=0, so it localises the 2-D search "
            "WITHOUT removing a single zero of the 1-D load (measured: "
            "sign changes of L_D inside the admissible arcs = sign changes "
            "of L_D over the whole period, at all 30 face witnesses)."),
        "route_B_measurement": (
            "sign-law-pruned branch-and-bound (face_bb_signlaw.py, "
            "best-first by area) at F4 witness (beta=0.30, y=-0.50): live "
            "undecided area 19.61 -> 0.076 (0.39%) in 30000 steps WITH the "
            "prefilter, versus 1.74 (8.8%) in 40000 steps with the same "
            "ordering and the prefilter OFF, versus 18.68 (95%) after "
            "400000 LIFO steps in the pre-sign-law face_bb.py.  Order-only "
            "exclusion measured over the 28+2 witnesses: 50-91% of the "
            "(theta1, D) torus for the positive faces F2-F7, 5-49% for the "
            "mixed faces F1/collar."),
        "blocking_lemma": (
            "branch cap = all-line crossing bound of the load curve: for "
            "every teacher, every D, every (alpha,gamma) != 0, "
            "alpha*(s0*W(t,D)+s1*W(t-beta,D)) + "
            "gamma*(s0*h(t)+s1*h(t-beta)) has <= 4 sign changes per "
            "period (measured tight; the student load is one line); plus "
            "the D-profile crossing bound.  Refuted routes with slack "
            "numbers recorded in branch_cap_scan.py (Sturm slack 6, "
            "per-arc disconjugacy slack 20, per-arc ratio monotonicity "
            "and winding bookkeeping refuted)"),
    }
    json.dump(out, open("count_match.json", "w"), indent=1)
    print("wrote count_match.json")


if __name__ == "__main__":
    main()
