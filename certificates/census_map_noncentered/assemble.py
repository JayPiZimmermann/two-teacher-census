"""Assemble certificate.json (+ CERTIFICATE.md via make_md.py) from the run
artifacts: validation.json, coincident_*.json, separated_*.json."""
import glob
import json
import os

import mpmath
from mpmath import iv

import noncentered as J
import certify_coincident as C

HERE = os.path.dirname(os.path.abspath(__file__))
J.set_prec(200)


def load(pat):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        with open(f) as fh:
            out.append((os.path.basename(f), json.load(fh)))
    return out


def nstr(x, n=25):
    return mpmath.nstr(x, n)


def main():
    with open(os.path.join(HERE, "validation.json")) as fh:
        val = json.load(fh)

    u, b1, b2 = C.certify_beta_star()
    PIv = J.PI_IV()

    coin = load("coincident_*.json")
    sep = load("separated_*.json")

    # ---- coincident summary -------------------------------------------
    beta_partition = []
    unresolved_beta = []
    lens_y_lo, lens_y_hi = None, None
    lens_boxes_total = 0
    lens_wide = 0
    lens_beta_measure = mpmath.mpf(0)
    ysum_check = True
    for name, d in coin:
        cnts = {}
        for b in d["boxes"]:
            k = (b["Wtau_total"], b["Wpot_total"])
            cnts[str(k)] = cnts.get(str(k), 0) + 1
        beta_partition.append({
            "interval": d["interval"],
            "beta_range": d["beta_range"],
            "certified_boxes": d["certified_boxes"],
            "measure_of_range": d["measure_of_range"],
            "measure_covered": d["measure_covered"],
            "measure_uncovered_inside_range": d["measure_uncovered_inside_range"],
            "root_counts_(Wtau,Wpot)_over_one_2pi_period": cnts,
            "seconds": d["seconds"],
        })
        for uu in d["uncovered_boxes"]:
            unresolved_beta.append({"interval": d["interval"], **uu})
        for b in d["boxes"]:
            for lr in b["lens_roots"]:
                lens_boxes_total += 1
                if lr["y"] is None or lr["y_width"] is None or lr["y_width"] > 1e-2:
                    lens_wide += 1
                    continue
                a, c = mpmath.mpf(lr["y"][0]), mpmath.mpf(lr["y"][1])
                lens_y_lo = a if lens_y_lo is None else min(lens_y_lo, a)
                lens_y_hi = c if lens_y_hi is None else max(lens_y_hi, c)
            if b["lens_roots"]:
                lens_beta_measure += mpmath.mpf(b["width"])

    # ---- separated summary --------------------------------------------
    sep_windows = []
    for name, ws in sep:
        for w in ws:
            ys = [x["y"] for x in w["branch_boxes"] if x.get("y")]
            ylo = min((mpmath.mpf(a) for a, b in ys), default=None)
            yhi = max((mpmath.mpf(b) for a, b in ys), default=None)
            sep_windows.append({
                "beta": w["beta"], "centre": w.get("centre"),
                "delta_collar": w["delta_collar"],
                "min_box_width": w["min_box_width"],
                "half_domain": w["half_domain"],
                "half_domain_area": w["half_domain_area"],
                "excluded_area": w["excluded_area"],
                "certified_branch_area": w["certified_branch_area"],
                "unresolved_area": w["unresolved_area"],
                "closed_fraction_of_half_domain":
                    w["closed_fraction_of_half_domain"],
                "n_branch_boxes": w["n_branch_boxes"],
                "n_unresolved": w["n_unresolved"],
                "n_evaluations": w["n_evaluations"],
                "seconds": w.get("seconds"),
                "y_hull_of_certified_branch_boxes":
                    None if ylo is None else [nstr(ylo, 18), nstr(yhi, 18)],
                "branch_boxes_sample": w["branch_boxes"][:40],
                "unresolved_sample": w["unresolved"][:20],
            })

    out = {
        "object": "certified enumeration of the boundary set of the NONCENTERED "
                  "(plain-ReLU) two-student census map",
        "model": {
            "source": "website/widgets.js kernelOf(\"noncentered\"), lines "
                      "231/232/233 (phiJ, HJ, dHJ), 238 (kernelOf), 294 "
                      "(slopeAtom), 295-297 (wronskian), 946-1005 "
                      "(separatedScan), 393-398 (ratioCoord)",
            "phiJ": "(pi-x) cos x + sin x, x = mod(t,2pi) folded by x>pi -> 2pi-x",
            "HJ": "|pi-x| sin x, x = mod(t,2pi)",
            "slopeAtomJ": "phiJ(t) - 2|sin t| = HJ'(t)",
            "kappa": "pi", "period": "2pi",
        },
        "precision": {
            "coincident_stratum_bits": C.PREC,
            "separated_stratum_bits": 90,
            "library": "mpmath %s interval arithmetic (directed rounding)"
                       % mpmath.__version__,
        },
        "validation": val,
        "beta_star": {
            "characterisation": "tan u = pi - u with u in (0, pi/2); "
                                "beta*_1 = 2u*, beta*_2 = 2pi - beta*_1",
            "u_star": [nstr(J.lo(u), 30), nstr(J.hi(u), 30)],
            "beta_star_1": [nstr(J.lo(b1), 30), nstr(J.hi(b1), 30)],
            "beta_star_2": [nstr(J.lo(b2), 30), nstr(J.hi(b2), 30)],
            "sum_minus_2pi": [nstr(x, 8) for x in
                              J.endpoints(b1 + b2 - 2 * PIv)],
        },
        "coincident_stratum": {
            "beta_partition": beta_partition,
            "unresolved_beta_boxes": unresolved_beta,
            "lens": {
                "certified_root_boxes": lens_boxes_total,
                "boxes_with_wide_or_undetermined_y": lens_wide,
                "beta_measure_carrying_the_lens": float(lens_beta_measure),
                "y_hull_of_tight_boxes":
                    None if lens_y_lo is None else [nstr(lens_y_lo, 20),
                                                    nstr(lens_y_hi, 20)],
            },
        },
        "separated_stratum": {
            "scope": "certified beta WINDOWS, not a covering of the beta axis",
            "windows": sep_windows,
        },
    }

    with open(os.path.join(HERE, "certificate.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("wrote certificate.json")
    return out


if __name__ == "__main__":
    main()
