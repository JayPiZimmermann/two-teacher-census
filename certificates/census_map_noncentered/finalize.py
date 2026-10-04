"""Build the scoped evidence artifact ``certificate.json``.

The artifact contains pointwise/scoped results and finite sample groupings.  It
does not certify the proposed global arrangement; ``CERTIFICATE.md`` is the
dated authoritative claim ledger.
"""
import collections
import glob
import json
import os
import subprocess

import mpmath
from mpmath import iv

import noncentered as J
import certify_coincident as C

HERE = os.path.dirname(os.path.abspath(__file__))
J.set_prec(200)


def jload(pat):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        with open(f) as fh:
            out.append((os.path.basename(f), json.load(fh)))
    return out


def flat(pat, key=None):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        with open(f) as fh:
            d = json.load(fh)
        out += d if isinstance(d, list) else [d]
    return out


def n(x, k=25):
    return mpmath.nstr(x, k)


# ---------------------------------------------------------------- inputs
with open(os.path.join(HERE, "validation.json")) as fh:
    VAL = json.load(fh)
u, B1, B2 = C.certify_beta_star()
PIv = J.PI_IV()
COIN = [(n_, d) for n_, d in jload("coincident_[0-9].json")]
SEP2D = flat("separated_*.json")
MAP = [p for p in flat("map_*.json") if "census" in p]
BAND = [p for p in flat("band_*.json") if "census" in p]
YW = [p for p in flat("ywall_*.json") if "census" in p]
CMP = None
if os.path.exists(os.path.join(HERE, "band_compare.json")):
    with open(os.path.join(HERE, "band_compare.json")) as fh:
        CMP = json.load(fh)

# ---------------------------------------------------------------- coincident
part = []
uncov = []
lens_lo = lens_hi = None
lens_boxes = 0
for name, d in COIN:
    cnt = collections.Counter()
    for b in d["boxes"]:
        cnt[(b["Wtau_total"], b["Wpot_total"])] += 1
    part.append({
        "interval": d["interval"], "beta_range": d["beta_range"],
        "certified_boxes": d["certified_boxes"],
        "measure_of_range": d["measure_of_range"],
        "measure_covered": d["measure_covered"],
        "measure_uncovered_inside_range": d["measure_uncovered_inside_range"],
        "counts_Wtau_Wpot": {str(k): v for k, v in cnt.items()},
        "seconds": d["seconds"],
    })
    for uu in d["uncovered_boxes"]:
        uncov.append({"interval": d["interval"], "beta": uu["beta"],
                      "width": uu["width"]})
    for b in d["boxes"]:
        for lr in b["lens_roots"]:
            lens_boxes += 1
            if lr["y"] is None or (lr["y_width"] or 1) > 1e-2:
                continue
            a, c = mpmath.mpf(lr["y"][0]), mpmath.mpf(lr["y"][1])
            lens_lo = a if lens_lo is None else min(lens_lo, a)
            lens_hi = c if lens_hi is None else max(lens_hi, c)

# ------------------------------------------------ sampled count-change pairs
sampled_y_count_change_brackets = {}
for b in sorted({q["beta"] for q in YW}):
    row = sorted([q for q in YW if q["beta"] == b], key=lambda q: -q["y"])
    lo_ = hi_ = None
    prev = None
    for q in row:
        if not q["certified"]:
            continue
        v = q["n_spur"]
        if prev is not None and v != prev[0]:
            hi_, lo_ = prev[1], q["y"]
        prev = (v, q["y"])
    if lo_ is not None:
        sampled_y_count_change_brackets[b] = {
            "y_between": [lo_, hi_],
            "n_spurious_above": 1,
            "n_spurious_below": 0,
            "scope": "adjacent certified samples; no intervening wall certified",
        }

# ----------------------------------------------- finite map-sample groupings
sample_points_by_census = collections.OrderedDict()
for p in MAP:
    if not p["certified"]:
        continue
    sample_points_by_census.setdefault(p["census"], []).append(
        (p["beta"], p["y"]))

OUT = {
    "object": "pointwise certified boundary/census data and a proposed "
              "arrangement of the NONCENTERED (plain-ReLU) two-student "
              "census map; global wall connectivity and face counts are "
              "not certified",
    "evidence_scope": "Certified entries are pointwise or scoped box results. "
                      "Grouping samples by proposed labels and interpolating "
                      "transition brackets do not certify global faces, walls, "
                      "continuation, or Euler counts.",
    "authoritative_ledger": {
        "path": "CERTIFICATE.md",
        "status_date": "2026-08-24",
    },
    "model": {
        "source": "website/widgets.js kernelOf(\"noncentered\") (line 238); "
                  "atoms 231/232/233, slopeAtom 294, wronskian 295-297, "
                  "Pot/Trq 240/241, separatedScan 946-1005, schurLabel 1202, "
                  "coincidenceTypeAt 1019-1027, ratioCoord 393-398",
        "phiJ": "(pi-x) cos x + sin x, x = mod(t,2pi) folded by x>pi -> 2pi-x",
        "HJ": "|pi-x| sin x, x = mod(t,2pi)",
        "slopeAtomJ": "phiJ(t) - 2|sin t| = HJ'(t)",
        "kappa": "pi", "period": "2pi",
        "kink_lattice_in_t": "{0, beta, pi, beta+pi} mod 2pi (four points; the "
                             "pair {pi, beta+pi} is the GHOST KINK with no "
                             "centered counterpart)",
    },
    "precision": {
        "coincident_bits": C.PREC,
        "per_teacher_census_bits": 160,
        "two_dimensional_sweep_bits": "140 (value) / 90 (sweep)",
        "library": "mpmath %s interval arithmetic, directed rounding"
                   % mpmath.__version__,
    },
    "validation": VAL,
    "reduction": {
        "statement": "Wtau = |pq| (sin d + d sinc(p) sinc(q)), "
                     "Wpot = |pq| (-sin d + d sinc(p) sinc(q)), "
                     "Wwgt = |pq| d sinc(p) sinc(q), with p = pi - (t mod 2pi) "
                     "folded to [-pi,pi], q = pi - ((t-beta) mod 2pi), "
                     "d = q - p in beta + 2pi Z",
        "consequence": "Wtau + Wpot = 2 Wwgt identically",
        "kinks": "absorbed into the four piece boundaries p,q in {0, +-pi}",
    },
    "symmetry": {
        "statement": "Wtau(2pi-beta, 2pi-t) = -Wtau(beta,t) and likewise for "
                     "Wpot, Wwgt; the mass direction is negated, which is the "
                     "(-,-) fold, so y is INVARIANT",
        "consequence": "fundamental domain beta in [0, pi]",
        "measured_on_shipped_classifier": "identical at 1536 of 1536 grid "
                                          "teachers (48 x 32 over the full map)",
    },
    "beta_star": {
        "characterisation": "tan u = pi - u, u in (0,pi/2); beta*_1 = 2u*, "
                            "beta*_2 = 2pi - beta*_1",
        "u_star": [n(J.lo(u), 30), n(J.hi(u), 30)],
        "beta_star_1": [n(J.lo(B1), 30), n(J.hi(B1), 30)],
        "beta_star_2": [n(J.lo(B2), 30), n(J.hi(B2), 30)],
        "sum_minus_2pi": [n(x, 8) for x in J.endpoints(B1 + B2 - 2 * PIv)],
    },
    "coincident_stratum": {
        "root_counts": {
            "Wpot": "constant 2 on (0,2pi)\\{pi}: t = pi and t = beta+pi only",
            "Wtau": "2 outside [beta*_1, beta*_2], 4 inside, 3 (a double root) "
                    "at the two endpoints, 2 at beta = pi",
            "Wwgt": "4 (t in {0, beta, pi, beta+pi}); determined by the "
                    "identity Wtau + Wpot = 2 Wwgt",
        },
        "beta_partition": part,
        "uncovered_beta_boxes": uncov,
        "lens": {
            "certified_root_boxes": lens_boxes,
            "y_hull_of_tight_boxes": None if lens_lo is None
            else [n(lens_lo, 20), n(lens_hi, 20)],
            "y_symmetry": "y_upper + y_lower = -1 exactly (t -> beta - t)",
        },
    },
    "separated_stratum": {
        "two_dimensional_mass_free_sweep": SEP2D,
        "per_teacher_certified": {
            "method": "locate (dense float grid, geometric in D) then certify "
                      "(2x2 Krawczyk in (s,D) = (th1, th0-th1))",
            "well_conditioned_direction": "weights from the RADIAL rows, "
                                          "massDet = pi^2 - phiJ(D)^2, which "
                                          "vanishes only at D = 0 (mod 2pi); "
                                          "widgets.js divides by h(D) instead "
                                          "and must guard |h(D)| > 5e-3, which "
                                          "blinds it to the antipodal gap",
            "type_test": "widgets.js schurLabel in interval arithmetic",
        },
    },
    "map_scan": {"points": MAP, "sampled_points_by_census":
                 {k: v for k, v in sample_points_by_census.items()}},
    "band_scan": {"points": BAND, "comparison": CMP},
    "sampled_y_count_change_brackets": sampled_y_count_change_brackets,
}

with open(os.path.join(HERE, "certificate.json"), "w") as fh:
    json.dump(OUT, fh, indent=1)
print("wrote certificate.json  (%d map points, %d band points, "
      "%d sampled y count-change brackets)"
      % (len(MAP), len(BAND), len(YW)))
