# The certified centered census arrangement and visual face pairing

This directory certifies the open-face arrangement and compares its signatures
with the source classifier and finite-grid renderer.  It does not infer exact
connected components from the renderer's flood fill.

It answers it in three links, and only the middle one is a computation.

| link | what it says | where |
|---|---|---|
| **containment** | the hypothesis-free classification of every centered two-student critical configuration and its exact variational type; the older positive- and mixed-sign stratum theorems remain ingredients and cross-checks | Lean; `LeanFormalization.Main.centered_two_student_full_classification`, with the ingredients and exact hypotheses in `ARRANGEMENT.md` section 7 |
| **enumeration** | the nontrivial interior zero sets give the two curves below, with every such root enclosed at 200 bits; the weight determinant also supplies fixed horizontal relabelling strata, and degenerate boundary columns are handled separately | `cert/`, `cert/CERTIFICATE.md` |
| **arrangement** | the two curves together with those fixed strata and boundary columns cut the map into exactly 10 faces, Euler-checked, each with a certified interior point and a certified census | `ARRANGEMENT.md`, `arrangement.json` |

## The headline

| certified or rendered object | count |
|---|---|
| faces of the arrangement | **10** |
| distinct census values on those faces | **4** |
| explicitly labelled, **non-topological visual face pairs** | **5** |
| finite-grid flood-fill groups drawn by the shipped map (`census-map.js` at commit `2201dd0`) | **5** |

Each visual pair has equal certified face signatures and matching listed local
exact-column witnesses.  This is a presentation relation, not a connectivity,
adjacency, or path-exclusion theorem.  Separately, each of the five finite-grid
renderer representatives reproduces the certified census of its face, string
for string (`face_widget.json`).  The renderer did not always have this output:
the asset of commit `72fb166` drew eight flood-fill groups — four tolerance
artifacts and a coalesced equal-signature grouping across a pinch no grid can
resolve.  `ARRANGEMENT.md`
sections 0 and 4 tell that story in full; the website fixes landed in
commits `c975b0e` (scale-relative selector) and `f430fab` (the certified
curve as a separator).

The exact orthogonal column matters even though a grid never samples it.
At `beta = pi/2`, both torque atoms share the roots `t = 0, pi/2`; the old
open-period scan omitted `t = 0` and falsely displayed a fit-only column.
The corrected classifier inserts both roots algebraically.  With
`c = (2/pi) atan(2/pi)`, it has one/two/one traps on each of
`(-1,0)` and `(0,1)`, with transition points at
`c-1,-c,c,1-c`.  Five selected column intervals provide local witnesses for
the labelled visual face pairs `F1|F4`, `F2|F3`, `F5|F8`, `F6|F9`,
`F7|F10`.  `face_widget.json` records the closed-form column cases, both
projective seams, source hashes, exact pair labels, and their agreement with
`website/widgets.js`; it makes no connected-component claim.

## The two curves

Both come from eliminating the teacher masses, which is possible because each
vanishing condition is two homogeneous linear equations in the mass pair
`(s0, s1)`.  What survives is a determinant in the gap and the angle alone.

* the **torque lens**, `torqueAtomWronskian = 0`, over
  `beta in [1.42092547555..., 1.72066717803...]` (`beta1*` solves
  `tan(beta/2) = pi/2 - beta/2` and `beta2* = pi - beta1*`);
* the **potential curve**, `potentialAtomWronskian = 0`, spanning the whole
  of `0 < beta < pi` with certified `y`-range `(c, 1-c)`,
  `c = (2/pi) arctan(2/pi) = 0.36090707...`.

The third determinant is not independent as a function:
`Wtau + Wpot = 2 Wwgt`, proved in Lean as
`torqueAtomWronskian_add_potentialAtomWronskian`
(`Planar/GeneralTeachers/PotentialBoundary.lean`).  This identity reduces the
number of independent determinant functions, but it does not identify their
zero loci.  The two nontrivial interior curves, the fixed horizontal
`Wwgt = 0` relabelling strata, and the separately treated degenerate columns
together supply the arrangement.

## The visual grouping the map used to coalesce

The potential curve has two distinct branches across `0 < beta < pi`.  The
outer face bands `0 < y < y_lower` and `y_upper < y < 1` carry the *same*
census, while their sampled separation pinches below one cell near the
degenerate boundary columns.  An unlabelled flood fill therefore coalesces
their samples; this is a finite-grid effect, not evidence about exact
connected components.

Refinement cannot draw this.  As `beta -> 0` the lens width tends to
`0.28113816169 * beta` (50-digit computation; the earlier quoted value
`0.2811381619` was wrong in its last digits), while the first column of a
depth-`d` grid sits at `beta = pi/(48 * 2^d)` and its cells are
`1/(6 * 2^d)` tall — both scale like `2^-d`, so the pinch is `0.110403`
cells wide at EVERY depth (measured constant to six digits over depths
4-16).  `website/precompute_census_map.js` therefore injects the certified
curve as a separator; `ARRANGEMENT.md` section 4 has the mechanism and the
guard rails (the torque lens must NOT be used as one).

## Containment

`ARRANGEMENT.md` section 7 states the Lean containment with exact hypotheses.
The public endpoint
`LeanFormalization.Main.centered_two_student_full_classification` has no
hypotheses on teacher or student angles or masses: criticality is equivalent to
the exhaustive labeled classification with its exact variational type.  The
older positive-mass theorem says that the local minima are exactly the global
minima plus the opposite-split coincidence traps at teacher-potential maxima
(`centered_twoTeacher_signed_localMin_iff_global_or_localMaxOuterSplit`).  For
mixed-sign teachers, the separated live stratum reduces to the exact fit
(`centered_mixedSign_distinct_critical_isExactFit`), quarter-gap configurations
do not exist (`centered_mixedSign_no_quarter_critical_point`), and the
equal-mass quarter family is never a local minimum
(`genQuarterBranch_eqTeacherMass_not_isLocalMin`).  The hypothesis-free public
endpoint also covers mixed-sign coincident configurations and zero-mass rows;
there is no residual row-typing obligation outside Lean.

## Reproducing

Requires `mpmath` (and `node` for the widget comparisons).  From this
directory:

```
node cert/ref_widgets.js           # reference samples from the literal widgets.js
python3 cert/validate.py           # kernel validation (writes cert/validation.json)
python3 cert/certify.py            # the boundary certificate  (long)
python3 cert/make_md.py            # regenerates cert/CERTIFICATE.md
python3 build_arrangement.py       # faces, Euler check, per-face census
node check_map.js                  # shipped classifier + asset vs the faces
python3 assemble.py                # regenerates arrangement.json
```

`cert/certificate_full.json.gz` (105 MB, every box of the exclusion search) is
**not** in version control; `cert/certify.py` regenerates it.  The compacted
`cert/certificate.json` kept here carries every enclosure and verdict.
`ref_samples.txt` (25 000 reference rows) is likewise regenerated by
`ref_widgets.js` rather than committed.

The shipped certificate is the **depth-16** run (`MAX_DEPTH_Y=16`), finished
2026-08-15.  `cert/PROVENANCE.json` records its run parameters, its console
log `cert/run_d16_shipped.log`, and the SHA-256 of both the compact and the
full artifact; it exists because that run predates the change that stamps the
parameters into every emitted certificate.

All enclosures use `mpmath.iv` at 200 bits with directed rounding
(`mp.prec = 220`); second-variation tests run at 60 decimal digits.
`float64` appears only as a locator for a high-precision solve and never
carries a verdict.
