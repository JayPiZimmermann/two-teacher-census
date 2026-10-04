"""Historical first-generation CERTIFICATE.md assembler.

The committed per-run shards needed by this script are no longer present. Its
arrangement section is a proposal assembled from finite witness/scanning data,
not a global wall or face-count certificate. The maintained CERTIFICATE.md
is the authoritative claim ledger (status 2026-08-24) and states the narrower
current evidence boundary.
"""
import collections
import glob
import json
import os

import mpmath

HERE = os.path.dirname(os.path.abspath(__file__))


def flat(pat):
    out = []
    for f in sorted(glob.glob(os.path.join(HERE, pat))):
        with open(f) as fh:
            d = json.load(fh)
        out += d if isinstance(d, list) else [d]
    return out


def jload(name):
    p = os.path.join(HERE, name)
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


# ------------------------------------------------------------ 3.4 coverage
rows = []
tot_unc = 0.0
for i, nm in ((0, "(0, beta*_1)"), (1, "(beta*_1, pi)"),
              (2, "(pi, beta*_2)"), (3, "(beta*_2, 2pi)")):
    d = jload("coincident_%d.json" % i)
    if d is None:
        rows.append("| `%s` | *not completed in this run* | | | |" % nm)
        continue
    cnt = collections.Counter((b["Wtau_total"], b["Wpot_total"])
                              for b in d["boxes"])
    tot_unc += d["measure_uncovered_inside_range"]
    rows.append("| `%s` | %d | %.15f | %.4g | %s |"
                % (nm, d["certified_boxes"], d["measure_covered"],
                   d["measure_uncovered_inside_range"],
                   ", ".join("%s on %d boxes" % (k, v)
                             for k, v in sorted(cnt.items()))))
cov = ["| beta interval | certified boxes | measure covered | uncovered inside | certified `(Wtau, Wpot)` counts |",
       "|---|---|---|---|---|"] + rows
cov.append("")
cov.append("Total uncovered beta measure inside the four intervals: "
           "**%.4g**.  The uncovered part consists only of collars around the "
           "critical betas `beta*_1`, `pi`, `beta*_2` (and, for the outer two "
           "intervals, around `beta = 0` and `beta = 2pi`), where `sin beta` "
           "or the piece structure degenerates; the critical betas themselves "
           "are treated exactly in section 3.2 and 3.3." % tot_unc)
sp = jload("coincident_spotcheck.json")
if sp:
    cov.append("")
    cov.append("The two outer intervals `(0, beta*_1)` and `(beta*_2, 2pi)` "
               "carry no lens and their counts are `(Wtau, Wpot) = (2, 2)` "
               "throughout.  A certified spot check on both (`certify_beta_box` "
               "on a shrinking box about each centre, accepting the widest box "
               "that certifies):")
    cov.append("")
    cov.append("| beta centre | certified half-width | `(Wtau, Wpot)` |")
    cov.append("|---|---|---|")
    for z in sp:
        if z["ok"]:
            cov.append("| %s | %.3g | (%d, %d) |" % (z["beta_centre"],
                                                     z["halfwidth"], z["Wtau"],
                                                     z["Wpot"]))
        else:
            cov.append("| %s | *not certified at any width down to* "
                       "`dist/2^19` | -- |" % z["beta_centre"])
    nok = sum(1 for z in sp if z["ok"])
    cov.append("")
    cov.append("**%d of %d** centres certify, all with `(2, 2)`.  The one that "
               "does not is `beta = 0.01`, within `0.01` of the degenerate "
               "column `beta = 0` where all three determinants vanish "
               "identically." % (nok, len(sp)))
COV = "\n".join(cov)

# ------------------------------------------------------------ 4.2 sweep
SEP2D = """The two-dimensional mass-free sweep was run as an interval branch and bound in
the `(s, D) = (th1, th0 - th1)` half domain `s in [0, 2pi]`, `D in [delta,
pi - delta]` (the whole stratum: of `D` and `2pi - D` exactly one is `<= pi`,
and `Dsep` is antisymmetric), with `Psi = Dsep/h(D)^4` for exclusion and an
interval Newton test in whichever coordinate has a nonvanishing partial for
existence and uniqueness of a branch.

**Measured conditioning** (thin `beta = 1.22`, uniform boxes over the half
domain with `delta = 0.15`, exclusion by `Psi` with the mean-value refinement
and the algebraically simplified gradient):

| box side `w` | boxes | excluded by `Psi` |
|---|---|---|
| 0.100 | 1827 | 25.3% |
| 0.050 | 7182 | 50.7% |
| 0.025 | 28728 | 67.0% |

On a `w = 0.025` box the enclosure of `Dsep` is a **median 10.7 times** wider
than its true range on that box, and the enclosure of `dDsep/dth0` a median
**52 times** wider; a float sign map at that resolution shows only **5.7%** of
cells genuinely contain the curve, so the exclusion is a factor of about five
short of what the geometry allows.  Adding the `beta`-width of a covering box
to the enclosure on the same footing as the two angle widths makes a covering
of the `beta` axis at a useful angular resolution unreachable at this cost.

**What this means for the certificate.**  The mass-free surface
`{separatedDet = 0}` in `(beta, th0, th1)` is NOT certified here as a complete
enumeration.  What IS certified is the per-teacher enumeration of section 4.3,
which answers the question the census actually asks.  The honest statement is:
*the two-dimensional exclusion closes about two thirds of the torus at the
finest resolution tried, and the remaining third is reported as uncovered, not
assumed empty.*"""

# ------------------------------------------------------------ ywall
YW = [p for p in flat("ywall_*.json") if "census" in p]
lines = ["| `beta` | sampled count change lies between | `n` spurious separated above | below |",
         "|---|---|---|---|"]
for b in sorted({q["beta"] for q in YW}):
    row = sorted([q for q in YW if q["beta"] == b and q["certified"]],
                 key=lambda q: -q["y"])
    prev = None
    for q in row:
        if prev is not None and q["n_spur"] != prev["n_spur"]:
            lines.append("| %.2f | `(%.3f, %.3f)` | %d | %d |"
                         % (b, q["y"], prev["y"], prev["n_spur"], q["n_spur"]))
            break
        prev = q
lines.append("")
lines.append("The certified endpoint samples bracket a change from 3 to 1 "
             "separated families. This is consistent with a saddle-node "
             "candidate, but the scan does not certify fold transversality, "
             "an intervening wall, or continuation in beta.")
YWALL = "\n".join(lines)

# ------------------------------------------------------------ arrangement
FA = jload("faces_all.json") or []
face_census = {}
for z in FA:
    face_census.setdefault(z["face"], []).append(z)


def census_of_face(fid):
    zs = face_census.get(fid, [])
    if not zs:
        return "?"
    c = collections.Counter(z["census"] for z in zs)
    return c.most_common(1)[0][0]


ARR = """**Domain.**  `beta in [0, pi]` is a fundamental domain: the census is invariant
under `beta -> 2pi - beta` (section 2.4 -- proved from the kernel parities, and
measured identical at 1536 of 1536 sampled teachers).  `y` is a CIRCLE of
circumference 2, because `y = -1` and `y = +1` are the same line `s0 = 0` (the
`(-,-)` fold).  So the fundamental domain is a **closed cylinder / annulus**,
whose two boundary circles are the columns `beta = 0` and `beta = pi`.  Both
are genuine strata of the arrangement, not merely cuts.

**A second exact symmetry.**  Swapping the two teacher atoms is a rotation
composed with `beta -> 2pi - beta`, so it acts on the map by `y -> -1-y`
(equivalently `psi -> pi/2 - psi`) at fixed `beta`.  Together with the
`beta`-mirror this makes the whole map symmetric about `y = -1/2` in the
mixed-sign sector and about `y = +1/2` in the same-sign sector -- which is
exactly the pairing `F2/F6` and `F3/F7` below.

**Proposed strata and sampled transitions.** The analytic fixed strata and
regional torque-lens certificates are combined in the proposed drawing with
sampled candidate transitions 1D-12 and 1D-13. The `7 x 10` scan down to
`y=0.005` records the same spurious-family count at its certified sample
points; it does not prove that no wall lies between them.

**Proposed arrangement bookkeeping on the fundamental-domain cylinder.**

```
   vertices  V = 11
   edges     E = 18
   faces     F =  7
   V - E + F = 0 = chi(annulus)
```

If all proposed connections are accepted, the eleven vertices are `0D-1`
(the lens double root `(beta*_1, -1/2)`),
`0D-3`, `0D-4` (the two points where the lens terminates on `beta = pi`),
`0D-5`, `0D-6` (`beta = 0` meets `y = 0` and `y = +-1`), `0D-9`..`0D-12` (the
two separated walls meeting each boundary column) and the two crossings of the
lens with the separated walls.  The eighteen edges are: 2 + 2 for the two lens
branches (each cut once by a separated wall), 1 for `y = 0`, 1 for `y = +-1`,
2 + 2 for the two separated walls (each cut once by the lens), and 4 + 4 for
the two boundary columns (each a `y`-circle cut at four points).

**Under the same proposal on the full map** `beta in [0, 2pi)` with `beta = 0` identified with
`beta = 2pi`, the domain is a TORUS and the two mirror copies are glued along
the genuine stratum `beta = pi`:

```
   V = 14 ,  E = 28 ,  F = 14 ,   V - E + F = 0 = chi(torus)
```

**Seven proposed regions, each with certified witness points.** `n` is the
pointwise certified number of torque roots. Agreement among a region's listed
witnesses does not prove constancy throughout that region.

| face | region | interior point | `n` | certified census | shipped agrees |
|---|---|---|---|---|---|
"""
FACE_META = [
    ("F1", "`0 < beta < pi`, `0 < y < 1` (same-sign sector)", "(0.75, 0.40)"),
    ("F2", "`w_a(beta) < y < 0`, outside the lens", "(0.75, -0.02)"),
    ("F3", "`w_a(beta) < y < 0`, inside the lens", "(3.10, -0.06)"),
    ("F4", "`-1-w_a < y < w_a`, outside the lens", "(0.75, -0.40)"),
    ("F5", "`-1-w_a < y < w_a`, inside the lens", "(2.60, -0.50)"),
    ("F6", "`-1 < y < -1-w_a`, outside the lens", "(0.75, -0.98)"),
    ("F7", "`-1 < y < -1-w_a`, inside the lens", "(3.10, -0.94)"),
]
NROOT = {"F1": 2, "F2": 2, "F3": 4, "F4": 2, "F5": 4, "F6": 2, "F7": 4}
AGREE = {"F1": "yes", "F2": "yes", "F3": "**no** (band, section 8)",
         "F4": "yes", "F5": "yes", "F6": "yes",
         "F7": "**no** (band, section 8)"}
CENS = {
    "F1": "`coincident:trap@positive | coincident:trap@positive | fit:global | separate:trap@mixed | separate:trap@mixed`",
    "F2": "`coincident:trap@mixed | fit:global | separate:trap@positive`",
    "F3": "`coincident:trap@mixed | coincident:trap@mixed | fit:global | separate:trap@positive`",
    "F4": "`coincident:trap@mixed | fit:global`",
    "F5": "`coincident:trap@mixed | coincident:trap@mixed | fit:global`",
    "F6": "`coincident:trap@mixed | fit:global | separate:trap@positive`",
    "F7": "`coincident:trap@mixed | coincident:trap@mixed | fit:global | separate:trap@positive`",
}
for fid, reg, rep in FACE_META:
    ARR += "| %s | %s | `%s` | %d | %s | %s |\n" % (
        fid, reg, rep, NROOT[fid], CENS[fid].replace("|", "\\|"),
        AGREE[fid])

ARR += """
`w_a(beta)` denotes the interpolated candidate transition used by the proposed
drawing, and `beta_c` its proposed crossing with the torque lens's upper
branch. The sampled crossing is bracketed only by **`(2.60,3.05)`**: at `beta=2.60`
the lens spans `y in (-0.72, -0.28)` and does not reach the band
`w_a = -0.10 < y < 0`, so `F3` is empty there, while at `beta = 3.05` the point
`(3.05, -0.09)` is certified to be inside the lens (4 torque roots) and inside
the band (`separate:trap@positive` present).  `F3` and `F7` therefore exist and
have certified witness points, but no complete global region or beta extent is
certified.

**The three notions of "piece", separately.**

| notion | count (fundamental domain) | count (full map) |
|---|---|---|
| faces in the proposed arrangement | 7 | 14 |
| distinct census values on those faces | **5** | **5** |
| level-set components predicted by that proposal | 7 | 14 |
| pieces reported by the shipped mosaic (depth 8) | -- | **30 big + 1868 fragments** |

Under the proposal, the level-set components equal the face count because the only two faces that
share a census, `F2/F6` and `F3/F7`, are separated by the whole middle band
`F4/F5`: they are the `y -> -1-y` mirror images of each other and never touch.

**What the shipped mosaic's 30 + 1868 is made of.**  The 1868 sub-threshold
fragments are the antipodal-band locator failure of section 8 (1780 of them
have `beta/pi in [0.97, 1.03]`, 21 more sit at `beta` within `0.05` of `0` or
`2pi`). Of the 30 big pieces, 14 are identified with the seven proposed
regions doubled by the `beta`-mirror; the remainder are the same locator failure at a larger
scale, plus the measure-zero strata `beta = 0`, `beta = pi`, `y = 0`,
`y = +-1` rendered with positive width."""

# ------------------------------------------------------------ band
CMP = jload("band_compare.json")
if CMP:
    causes = CMP["mismatch_causes"]
    BAND = """**The measurement.**  `12 x 12` grid over `beta in [3.05, 3.235]`,
`y in (-1, 1)`; at each point the CERTIFIED census of section 5 (certified
sign chart plus individually Krawczyk-certified separated families) against the
SHIPPED `censusSignature("noncentered", .)` loaded from the unmodified
`website/widgets.js` through the `precompute_census_map.js` DOM shim.

```
   points                                   %d
   certified                                %d
   agree with the shipped signature         %d
   disagree                                 %d
```

**Every disagreement is of one sign.**  The mismatch causes are

```
   %s
```

read as `separate_<certified>_vs_<shipped>`: in **all %d** cases the shipped
classifier reports FEWER separated traps than are certified to exist.  Never
more, and never a coincident-count disagreement inside this band's certified
part.

**The certified sample values across the band.** On the 135 certified points
the census takes three values, grouped by sampled `y` as follows:

```
   y > 0                          e = coincident:trap@positive x2 | fit:global
                                      | separate:trap@mixed x2
   y < 0, |y| and |y+1| > w_a     a = coincident:trap@mixed x2 | fit:global
   y < 0, inside the two bands    b = a + separate:trap@positive
```

Among this finite grid, every deviation from that sampled pattern occurs at the
nine points the method declines to certify, all in the two `beta`
columns nearest `pi` (`beta = 3.1348` and `beta = 3.1502`).  The shipped map
over the same grid takes five values in an irregular, speckled pattern.

**Verdict: hypothesis (b).**  The fog is a FAILURE OF THE SHIPPED
`separatedScan` LOCATOR, not genuine fine structure.  The mechanism is
identified and is not a tolerance:

* The separated families that carry the census near `beta = pi` sit at student
  gaps `D = th0 - th1` of about **0.02 to 0.07**.  `separatedScan` seeds Newton
  from a `110 x 110` grid, spacing `0.057` -- one grid cell -- and only from
  seeds whose residual is already below `0.35`.  Near-coincident pairs fail
  that pre-filter intermittently, so the family is found at some `(beta, y)`
  and missed at neighbouring ones.  That is exactly a cell-scale flicker.
* `separatedScan` divides by `h(D)` and guards `|h(D)| > 5e-3`, which makes it
  blind to the entire ANTIPODAL stratum `D = pi`.  The certified enumeration
  solves the weights from the radial rows instead, whose determinant
  `pi^2 - phiJ(D)^2` vanishes only at `D = 0`, and does find separated families
  at `D = pi` exactly (a saddle at `beta = 2.6, y = -1/2`, for instance).
* Independently, `findRoots` (900-point scan on `t`, with a `1e-5`
  root-merging tolerance) misses one of the four torque roots when `beta` is
  within about `3e-3` of `pi`, because two pairs of roots collide at `t = 0`
  and `t = pi` there.  Measured at `beta = 3.14108`, `y = -0.03548`: certified
  4 roots and 2 coincident traps, shipped 1 trap.

**Comparison with the proposed count.** The proposed arrangement has 7 regions
on the fundamental domain; the mosaic reports 30 big pieces plus 1868
fragments. The observed concentration of fragments diagnoses locator behavior,
but it does not by itself certify the proposed global region count.

**The collar.** Outside `|beta-pi|>0.0154` the certified samples in the band
show no additional census value; inside it the method declines at some `y`.
The sampled statement is:

> at each certified sampled `beta` in `[3.05,3.235]` with
> `|beta-pi|>0.0154`, and each certified sampled `y`, the census takes one of
> the three listed values. No statement between grid points follows.

The `0.0154` is the width of the two `beta`-columns the scan could not
certify, not a proof that structure exists inside them; a finer scan there was
started and is reported as unresolved in section 9.""" % (
        CMP["n_points"], sum(1 for p in CMP["points"] if p["certified"]),
        CMP["n_match"], CMP["n_mismatch"],
        "\n   ".join("%-22s %d" % (k, v) for k, v in sorted(causes.items())),
        CMP["n_mismatch"])
else:
    BAND = "*band comparison not available*"

# ------------------------------------------------------------ unresolved
unc = []
for i in range(4):
    d = jload("coincident_%d.json" % i)
    if d:
        for u in d["uncovered_boxes"]:
            unc.append((d["interval"], u["beta"], u["width"]))
UN = ["### 9.1 beta collars where a coincident root COUNT is not certified", "",
      "| interval | box | width |", "|---|---|---|"]
for nm, b, w in unc:
    UN.append("| `%s` | `[%s, %s]` | %.3e |" % (nm, b[0], b[1], w))
UN.append("")
UN.append("Total: **%.4g** in beta measure.  Each collar sits around one of "
          "`0`, `beta*_1`, `pi`, `beta*_2`, `2pi`, all of which are treated "
          "exactly in sections 3.2-3.3." % sum(w for _, _, w in unc))
UN.append("")
UN.append("""### 9.2 the two-dimensional mass-free separated sweep

Not closed.  At the finest resolution tried the exclusion covers about **two
thirds** of the `(th0, th1)` half domain for a thin `beta`; the remaining third
is reported, not assumed empty.  See section 4.2 for the measured conditioning.
The `h(D) = 0` collar (`D = 0 (mod pi)`) is excluded by modelling and its area
is `2 * 2pi * delta` for the collar half-width `delta` used.

### 9.3 points at which the per-teacher census is not certified

The certified census declines when a determinant sign at a torque root
straddles zero, when the torque-root isolation leaves an undecided box, or when
a located separated candidate fails the Krawczyk test at every box radius down
to `1e-3 / 2^40`.  In the band scan that is **9 of 144** points, all in the two
`beta`-columns nearest `pi`.  In the `y > 0` wall scan it is the points with
`y <= 0.02` at `beta <= 1.0` and `y = 0.002` everywhere -- i.e. a collar around
the degenerate stratum `y = 0`, where the teacher has one effective atom.
In the fundamental-domain map scan (`map_scan.py`, stopped after 137 of 336
points to free compute for the band investigation) the uncertified points are
again exactly those with `|y|` or `|y| - 1` below about `0.05` at small `beta`.

### 9.4 completeness of the separated locator

The separated half of every certified census rests on a float locator
(`260 x 240` seeds, `D` sampled geometrically from `1e-4`).  Every family it
returns is certified individually; a family it never seeds is not excluded.
This is the single residual obligation of the face censuses, and the
two-dimensional sweep of section 4.2 is the instrument that would discharge it.

### 9.5 the same-sign sector near beta = 0

At `beta = 0.05` the certified census at `y = 0.02` and `y = 0.98` is
`coincident:trap@positive x2 | fit:global | separate:trap@mixed` -- one
separated trap, not two -- while at the sampled `beta>=0.15` points the count
is 2 down to `y=0.005`. This is evidence for a possible additional transition
below about `beta=0.1`; its existence, extent, and effect on any global face
count are unresolved.""")
UNRES = "\n".join(UN)

# ------------------------------------------------------------ write
p = os.path.join(HERE, "CERTIFICATE.md")
with open(p) as fh:
    md = fh.read()
for tag, txt in (("COINCIDENT_COVERAGE", COV), ("SEPARATED_2D", SEP2D),
                 ("YWALL", YWALL), ("ARRANGEMENT", ARR), ("BAND", BAND),
                 ("UNRESOLVED", UNRES)):
    md = md.replace("<!--%s-->" % tag, txt)
# The first-generation headline substitution is retired. The maintained
# CERTIFICATE.md carries an explicit scoped-evidence table and must not be
# overwritten by the old proposed-arrangement bookkeeping.
with open(p, "w") as fh:
    fh.write(md)
print("filled CERTIFICATE.md")
