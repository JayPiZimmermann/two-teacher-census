"""
Build the planar arrangement of the certified census boundary in the map
rectangle (beta, y), and the census on every face.

DOMAIN.  beta in [0, pi]; y in R/2Z (the mass-ratio coordinate is a CIRCLE:
y = -1 and y = +1 are the same line s0 = 0).  So the map is a closed CYLINDER,
not a rectangle, and beta = 0 / beta = pi are its two boundary circles.  (They
are also the same stratum of the kernel, since everything is pi-periodic in
beta; identifying them turns the cylinder into a torus without changing the
ten open faces -- see ARRANGEMENT.md.)

CURVES.  Certified in cert/CERTIFICATE.md:
  torque lens      image of {Wtau = 0}, beta in [beta1*, beta2*], y in (c-1,-c)
  potential curve  image of {Wpot = 0}, beta in (0, pi),          y in (c, 1-c)
  y = 0            image of {Wwgt = 0} at t = 0        (s1 = 0)
  y = +-1          image of {Wwgt = 0} at t = beta     (s0 = 0)
  beta = 0, beta = pi/2, beta = pi   degenerate columns

Face membership of every reported interior point is certified twice over:
  (a) by the certified y-enclosures of the curve branches at that beta, and
  (b) by the certified SIGN CHART at the point (root count + determinant signs),
      which is a complete invariant of the face.
"""
import json
import math
import mpmath
import cert_core as K

C = K.C
iv = K.iv

BETA1 = "1.420925475551033713494856535"     # beta*_1  (certificate section 6)
BETA2 = "1.720667178038759524968"           # beta*_2 = b*
CVAL = "0.360907073228108373304919068364"   # c = (2/pi) arctan(2/pi)


# ---------------------------------------------------------------------------
# the ten faces, with an interior point each (exact dyadic float64 values, so
# the interval inputs are degenerate and every enclosure below is certified)
# ---------------------------------------------------------------------------
FACES = [
    dict(id="F1", name="lower, left of beta=pi/2, outside the torque lens",
         region="0 < beta < pi/2,  -1 < y < 0,  outside the torque lens",
         rep=(0.75, -0.25),
         extra=[(0.4, -0.5), (1.5, -0.9), (0.25, -0.1), (1.5, -0.1)]),
    dict(id="F2", name="lower, left of beta=pi/2, INSIDE the torque lens",
         region="beta1* < beta < pi/2,  inside the torque lens",
         rep=(1.5, -0.5), extra=[(1.55, -0.45), (1.45, -0.5)]),
    dict(id="F3", name="lower, right of beta=pi/2, INSIDE the torque lens",
         region="pi/2 < beta < beta2*,  inside the torque lens",
         rep=(1.65, -0.5), extra=[(1.6, -0.55), (1.7, -0.5)]),
    dict(id="F4", name="lower, right of beta=pi/2, outside the torque lens",
         region="pi/2 < beta < pi,  -1 < y < 0,  outside the torque lens",
         rep=(2.5, -0.5), extra=[(1.65, -0.9), (2.9, -0.1), (1.8, -0.2)]),
    dict(id="F5", name="upper, left of beta=pi/2, BELOW the potential curve",
         region="0 < beta < pi/2,  0 < y < y_lower(beta)",
         rep=(0.4, 0.2), extra=[(1.5, 0.05), (0.05, 0.2), (1.2, 0.3)]),
    dict(id="F6", name="upper, left of beta=pi/2, INSIDE the potential lens",
         region="0 < beta < pi/2,  y_lower(beta) < y < y_upper(beta)",
         rep=(0.4, 0.5), extra=[(1.5, 0.5), (0.05, 0.5)]),
    dict(id="F7", name="upper, left of beta=pi/2, ABOVE the potential curve",
         region="0 < beta < pi/2,  y_upper(beta) < y < 1",
         rep=(0.4, 0.8), extra=[(1.5, 0.95), (0.05, 0.8), (1.2, 0.7)]),
    dict(id="F8", name="upper, right of beta=pi/2, BELOW the potential curve",
         region="pi/2 < beta < pi,  0 < y < y_lower(beta)",
         rep=(2.5, 0.2), extra=[(1.65, 0.05), (3.1, 0.2)]),
    dict(id="F9", name="upper, right of beta=pi/2, INSIDE the potential lens",
         region="pi/2 < beta < pi,  y_lower(beta) < y < y_upper(beta)",
         rep=(2.5, 0.5), extra=[(1.65, 0.5), (3.1, 0.5)]),
    dict(id="F10", name="upper, right of beta=pi/2, ABOVE the potential curve",
         region="pi/2 < beta < pi,  y_upper(beta) < y < 1",
         rep=(2.5, 0.8), extra=[(1.65, 0.95), (3.1, 0.8)]),
]

# These are presentation labels, not edges of an adjacency graph.  Each pair
# consists of two open faces with the same certified signature and names the
# selected exact beta=pi/2 source-classifier witnesses checked by check_map.js.  Neither
# this list nor that local witness establishes path connectivity or excludes a
# route elsewhere in the parameter space.
VISUAL_FACE_PAIR_SPECS = [
    dict(id="VP1", faces=["F1", "F4"],
         column_witness_cases=["lower_open_below_c_minus_1",
                               "lower_open_above_minus_c"],
         column_witness_intervals=["(-1,c-1)", "(-c,0)"]),
    dict(id="VP2", faces=["F2", "F3"],
         column_witness_cases=["lower_open_middle"],
         column_witness_intervals=["(c-1,-c)"]),
    dict(id="VP3", faces=["F5", "F8"],
         column_witness_cases=["upper_open_below_c"],
         column_witness_intervals=["(0,c)"]),
    dict(id="VP4", faces=["F6", "F9"],
         column_witness_cases=["upper_open_middle"],
         column_witness_intervals=["(c,1-c)"]),
    dict(id="VP5", faces=["F7", "F10"],
         column_witness_cases=["upper_open_above_1_minus_c"],
         column_witness_intervals=["(1-c,1)"]),
]

VERTICES = [
    ("0D-1", BETA1, "-0.5", "torque lens left endpoint (double root of Wtau)"),
    ("0D-2", BETA2, "-0.5", "torque lens right endpoint (double root of Wtau)"),
    ("0D-3", "pi/2", "-c", "torque lens upper branch meets beta = pi/2"),
    ("0D-4", "pi/2", "c-1", "torque lens lower branch meets beta = pi/2"),
    ("0D-5", "pi/2", "c", "potential lower branch meets beta = pi/2"),
    ("0D-6", "pi/2", "1-c", "potential upper branch meets beta = pi/2"),
    ("0D-7", "0", "0.5", "potential curve endpoint on beta = 0"),
    ("0D-8", "pi", "0.5", "potential curve endpoint on beta = pi"),
    ("0D-9", "0", "0", "beta = 0 meets y = 0"),
    ("0D-10", "0", "1", "beta = 0 meets the y = +-1 seam"),
    ("0D-11", "pi/2", "0", "beta = pi/2 meets y = 0"),
    ("0D-12", "pi/2", "1", "beta = pi/2 meets the y = +-1 seam"),
    ("0D-13", "pi", "0", "beta = pi meets y = 0"),
    ("0D-14", "pi", "1", "beta = pi meets the y = +-1 seam"),
]


def branch_ys(beta):
    """certified y-enclosures of the two curve branches present at beta."""
    out = {}
    for which in ("torque", "potential"):
        bs = K.curve_branches(beta, which)
        vals = []
        for b in bs:
            if b["y_lo"] is None:
                continue
            vals.append((mpmath.mpf(b["y_lo"]), mpmath.mpf(b["y_hi"])))
        vals.sort()
        out[which] = vals
    return out


def certify_point(beta, y):
    ch = K.sign_chart(beta, y)
    br = branch_ys(beta)
    Y = mpmath.mpf(y)
    pos = {}
    for which, vals in br.items():
        if not vals:
            pos[which] = "absent"
        elif len(vals) == 2:
            (lo1, hi1), (lo2, hi2) = vals
            if Y < lo1:
                pos[which] = "below"
            elif hi1 < Y < lo2:
                pos[which] = "inside"
            elif Y > hi2:
                pos[which] = "above"
            else:
                pos[which] = "UNDECIDED"
        else:
            pos[which] = "UNEXPECTED(%d)" % len(vals)
    return ch, {k: [[mpmath.nstr(a, 22), mpmath.nstr(b, 22)] for a, b in v]
                for k, v in br.items()}, pos


def main():
    out = {
        "object": "planar arrangement of the certified census boundary of the "
                  "CENTERED two-student census map",
        "domain": {
            "beta": "[0, pi]  (closed; beta = 0 and beta = pi are the same "
                    "stratum of the kernel, so the map is a cylinder and may be "
                    "closed to a torus without changing the face count)",
            "y": "R / 2Z  -- a CIRCLE; y = -1 and y = +1 are ONE line (s0 = 0)",
            "note": "psi = (y+1)pi/2, (s0,s1) = (sin psi, cos psi); the map's "
                    "representative always has s0 >= 0"},
        "topology": {
            "displayed_square_with_y_seam_identified": {
                "space": "closed annulus (equivalently, a cylinder with boundary)",
                "vertices": 14, "edges": 24, "faces": 10,
                "euler_characteristic": 0},
            "after_identifying_beta_0_with_beta_pi": {
                "space": "torus", "vertices": 11, "edges": 21, "faces": 10,
                "euler_characteristic": 0},
            "note": "beta = 0 and beta = pi are duplicate drawings of one "
                    "centered projective stratum; y = -1 and y = 1 are likewise "
                    "one mass-coordinate seam"},
        "precision": "mpmath.iv, 200 bits, directed rounding (mp.prec 220)",
        "strata_1d": [
            "beta = 0        (degenerate column, all three determinants vanish)",
            "beta = pi/2     (degenerate column, common zero set on the kink lattice)",
            "beta = pi       (degenerate column, = beta 0 by pi-periodicity)",
            "y = 0           (s1 = 0; image of the Wwgt root t = 0)",
            "y = +-1         (s0 = 0; image of the Wwgt root t = beta)",
            "torque lens     (image of Wtau = 0; closed curve, beta in [beta1*, beta2*])",
            "potential curve (image of Wpot = 0; two branches over beta in (0,pi))"],
        "vertices": [{"id": a, "beta": b, "y": c, "what": d}
                     for a, b, c, d in VERTICES],
        "faces": [],
    }

    for f in FACES:
        ch, br, pos = certify_point(*f["rep"])
        rec = {"id": f["id"], "name": f["name"], "region": f["region"],
               "interior_point": {"beta": repr(f["rep"][0]), "y": repr(f["rep"][1])},
               "certified": ch["certified"],
               "torque_root_count": ch["nroots"],
               "roots": ch["roots"],
               "curve_branch_y_enclosures_at_this_beta": br,
               "position_relative_to_curves": pos,
               "census": ch["census"]}
        wit = []
        for (bb, yy) in f["extra"]:
            c2, b2, p2 = certify_point(bb, yy)
            wit.append({"beta": repr(bb), "y": repr(yy),
                        "certified": c2["certified"],
                        "torque_root_count": c2["nroots"],
                        "position_relative_to_curves": p2,
                        "census": c2["census"]})
        rec["extra_witnesses"] = wit
        rec["census_constant_on_witnesses"] = all(
            w["census"] == rec["census"] for w in wit)
        out["faces"].append(rec)
        print("%-4s n=%d cert=%s pos=%s census=%s  (witnesses agree: %s)"
              % (f["id"], ch["nroots"], ch["certified"], pos, ch["census"],
                 rec["census_constant_on_witnesses"]), flush=True)

    cens = {}
    for r in out["faces"]:
        cens.setdefault(r["census"], []).append(r["id"])
    out["face_count"] = len(out["faces"])
    out["census_values"] = [{"census": k, "faces": v} for k, v in sorted(cens.items())]
    out["distinct_census_values"] = len(cens)
    face_by_id = {r["id"]: r for r in out["faces"]}
    pairs = []
    for spec in VISUAL_FACE_PAIR_SPECS:
        signatures = [face_by_id[face_id]["census"] for face_id in spec["faces"]]
        equal = len(set(signatures)) == 1
        if not equal:
            raise RuntimeError("visual face pair has unequal certified signatures: %r"
                               % spec)
        pairs.append(dict(spec, face_signatures=signatures,
                          shared_signature=signatures[0],
                          face_signature_agreement=equal))
    out["visual_face_pairs"] = {
        "count": len(pairs),
        "claim_scope": "explicit non-topological presentation pairs only; "
                       "no connectivity, adjacency, or path-exclusion claim",
        "connectivity_checked": False,
        "pairs": pairs,
    }
    with open("arrangement_faces.json", "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nfaces:", out["face_count"], " distinct censuses:", out["distinct_census_values"])
    for k, v in sorted(cens.items()):
        print("   %-70s %s" % (k, v))


if __name__ == "__main__":
    main()
