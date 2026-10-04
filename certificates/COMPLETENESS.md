# Completeness ledger for the census-map certificates

For each claim: the instrument, its precision, the exact domain on which the
claim is established, and — separately — what the instrument cannot see.
"Certified" = interval arithmetic / exact identity; "proved" = Lean theorem;
"measured" = float sampling.  Details and artifacts are in the two
directories; this file is the union of their honest scopes, kept in one
place so that no strengthening of one table silently overstates another.

Current status: 2026-08-24. Historical campaign details below are retained as
provenance; the two summary tables state the current result.

## Centered map (`census_map_centered/`)

| claim | status | domain | residual |
|---|---|---|---|
| root counts of the three mass-free determinants | certified (200-bit boxes + exact reduced-form inequalities C1–C5) | every `beta in (0, pi)` except the three values `beta*_1, pi/2, beta*_2` | the three values themselves (handled exactly, zero area; no interval method can certify AT a margin-degenerate point) |
| the two boundary curves (torque lens, potential curve): existence, uniqueness, `y`-enclosures | certified (200 bits) | lens over `[beta*_1, beta*_2]`, potential curve over `(0, pi)` | `y`-enclosures wider than `1e-3` on the boxes listed in cert §9.2 (certified, not tight); mass direction undetermined on 98 boxes adjacent to `beta = pi/2` |
| per-face census (10 faces, 35 witnesses = 10 face representatives + 25 further witnesses, counted from `arrangement.json` 2026-08-10) | certified sign charts | the listed witnesses; sign chart is a complete face invariant | census BETWEEN witnesses follows from the count constancy, itself certified per `beta`, not per point |
| quarter family on `y = -1/2` is never a minimum | proved (`genQuarterBranch_eqTeacherMass_not_isLocalMin`) | all `beta in (0, pi) \ {pi/2}` | interval corroboration stops at `0.002` from the columns; the Lean statement has no such gap |
| containment: classifier census = true census | proved by the hypothesis-free public endpoint `LeanFormalization.Main.centered_two_student_full_classification`, which packages the exhaustive labeled classification and `L_CVariationalType` | all two-student/two-teacher centered configurations; no sign, liveness, or direction hypotheses | numerical sign charts remain independent corroboration, not a residual typing obligation |
| shipped-map agreement (4 keys, 5 pieces, representatives) | measured (exact float inputs through the generator's DOM shim) | `census-map.js` and `widgets.js` at commit `2201dd0` | a statement about the shipped JavaScript at that commit, not about real numbers |

## Noncentered map (`census_map_noncentered/`)

| claim | status | domain | residual |
|---|---|---|---|
| coincident root counts (`Wtau`: 2/3/4; `Wpot`: 2) | certified (exact sympy identities + interval-certified inequalities, `analytic_counts.json`) | EVERY `beta in (0, 2pi) \ {pi}`; `beta = 0, pi` are exact degenerate cases | none — the analytic closure leaves no collar and no uncovered interval |
| `beta*` (lens endpoints) | certified enclosure, 28 digits | — | — |
| torque-lens branches, `y`-enclosures | certified (200-bit boxes) | `(beta*_1, pi)` and `(pi, beta*_2)`, 309 + 305 boxes | uncovered measure `2.0e-13` per interval (count itself closed analytically; only the branch ENCLOSURES rest on the boxes) |
| per-family existence and type | certified per family (2x2 Krawczyk at 160 bits; Schur type in interval arithmetic, equal to the tree's second-order test up to `massDet^2 > 0`) | every listed family with `D >= 1e-4` | outside the exhaustive teachers below, an unseeded family is not excluded; the float locator proposes candidates but is not completeness evidence |
| exact per-teacher enumerations | certified exhaustive branch-and-bound replay (`ivcert_*.bits.json`, `ivcert_check.py`) | 22 named teachers across the proposed labels F1–F7; 2,368,457 leaves, zero failed leaves and zero undecided boxes | pointwise only; no statement about intervening teachers or global cells |
| mass-free separated sweep | attempted, not closed | measured at thin `beta = 1.22`: exclusion 25/51/67% at box sides 0.1/0.05/0.025; enclosures 3.7× (value) and 17.9× (gradient) over true ranges; ~3.2% of finest cells genuinely contain the curve | interval dependency prevents this route from proving a global wall cover; it is no longer the evidence for completeness at the 22 exhaustive teachers |
| proposed-cell witness lists (7 labels, 28 witnesses) | complete coincident sign charts plus Krawczyk-certified separated families at each point; exact only at witnesses among the 22 exhaustive enumerations | `faces_fixed.json`; 26 original and 2 relocated witnesses, with minimum `\|beta-pi\| = 0.0416` | the interpolated fold arcs, their incidences, complement components, and a global seven-face arrangement are not certified; the labels organize sampled evidence |
| mixed-collar additional selected-row list | the coincident chart is complete and every located separated family is Krawczyk-certified at `(0.15, 0.004)` and its mirror; the separated complement is not covered, so this is not an exact census | nine sampled gaps `beta = 0.15, 0.30, 0.50, 0.75, 1.00, 1.50, 2.20, 2.60, 3.00`; no interval between them is claimed | only `beta = 0.15` has a two-sided locate-then-certify transition pair; the other eight rows are one-sided first-success heights. None proves a global wall, 1-D stratum, piece count, or collar-interior classification |
| proposed same-sign-sector fold traces (`w_a`) | complete coincident sign charts plus individually Krawczyk-certified returned-family lists on two sides | 8 sampled gaps in `[0.3, 3.13]`, with displayed windows of width `0.02`–`0.03` in `y` | these teachers lack exclusion covers for the separated complement; neither the pointwise windows nor the interpolating trace is a certified count boundary |
| mirror symmetry `beta -> 2pi - beta` | proved from kernel parities (an identity of the branch forms) | everywhere | — |
| shipped-classifier / locator-list agreement | measured, with per-returned-family interval checks | at `f430fab`: 94/146 accepted locator lists, one-sided failure, 1868 fragments; at `2201dd0`: 146/146 accepted lists and 49 components — **superseded by the Jacobian fix**, whose current comparison has agreement at **137/137 accepted locator lists** (156 of 158 points evaluated; `band_comparefix2.json`) | this is not an exhaustive census at the band points: outside the separately enumerated teachers an unlocated separated family is not excluded.  The other 19 evaluated points are precisely locator declines, all in the antipodal pocket.  The 2 unevaluated points (`beta = 3.129, 3.131` at `y = 0`) are one-teacher degeneracies (`s1 = 0`); the comparison makes no claim there |
| **the order-only prefilter CANNOT decide the Jacobian determinant** | certified counterexample (two Krawczyk-certified zeros at one teacher) | at `(beta, y) = (2.60, -0.480)` two separated critical configurations on the F5 continuation have the SAME sign-law order pattern `><><` and Jacobian determinants of OPPOSITE sign, `-86.2459` and `+40.1637`; the four compared offsets are `th0`, `th0-beta`, `-th1`, `beta-th1`, read off `generalJAngleEq0_unitTeacher0/1` and `generalJAngleEq1_unitTeacher0/1` (the first version of the test used `th1`, `th1-beta` and was corrected against the declarations; verdict unchanged).  Both tracked branches carry that one pattern at every certified teacher of the walk, `y = -0.500` through `-0.475`, with opposite determinant signs throughout | so region-level birth-freeness cannot be had from the sign law's order-only exclusion, and the count bootstrap stays circular by this route.  This is a NEGATIVE result about one method, not about the sign law (`SeparatedCount/SignLawJ.lean` is about `separatedNumJ` and is untouched) and not about an interval branch-and-bound on the determinant itself, which is the covering sweep and is priced at `~1e7` boxes per `0.002 x 0.002` cell |
| **continuation certificates on nine proposed-label rectangles** | certified (Krawczyk per step at 160 bits, Jacobian-determinant sign per step; emitted by `cont_cert.py`, replayed by `cont_check.py`, gated by `cont_gate.py`) | `cont_F1a`, `cont_F1b`, `cont_F2`, `cont_F3`, `cont_F4`, `cont_F5a`, `cont_F5b`, `cont_F6`, `cont_F7`; 2,589 grid teachers (2,572 distinct), 30,278 certified steps, zero undecided leaves | tracked-sheet persistence at grid teachers only; it neither excludes untracked zeros nor proves a global face arrangement, and it says nothing between grid teachers |
| **F5 point counts differ; continuation gives measured fold evidence** | exact at the two exhaustively enumerated teachers; certified only at each accepted continuation step; walk stopping and determinant-zero fits are prospecting | at `beta = 3.00`, the complete enumerations have three separated families at `y = -0.42` and one at `y = -0.40`. Guard-off walks last certify two tracked families at `y = -0.4980`, `-0.4580`, `-0.4240`, `-0.4160`, `-0.4140` for `beta = 2.30`, 2.60, 2.90, 3.00, 3.10, with opposite determinant signs at those accepted points. The instrument uses `0.002` steps, 160 bits, and disabled `maxjump`; it stops on `UNDECIDED_MINRAD`. The earlier guarded `beta = 3.00` value `-0.418` and the extrapolation `-0.4585` are withdrawn | the exact point counts prove only that the proposed F5 label cannot carry one constant count. The continuation trends are consistent with a fold, but a failed Krawczyk step is not a degeneracy certificate: no `separatedAngleJacDetJ = 0` point, fold curve, regional `hfold` verdict, or certified count-band subdivision follows |
| **all 22 per-teacher enumerations** | certified (sign-law-pruned branch and bound, packed and replayed by `ivcert_pack.py` / `ivcert_check.py`) | all teachers in the exact-enumeration row above; 2,368,457 leaves in total, including F1's 254,893-leaf witness | each result is complete at one teacher and does not imply completeness over its proposed cell label |
| **row 7a as a REGION statement; the naive mass-free exclusion REFUTED at the exact fit** | proved (`CountConstancyJ/SchurExclusionJ.lean`: 21 declarations, every theorem axiom-clean) | `censusStripJ seam delta dmax` with `0 < delta` and `dmax < pi`, over an arbitrary teacher set `S`; the exact-fit input is discharged outright whenever `delta` exceeds every teacher gap in `S` | the schema's `hfold` now follows from exactly two interval-testable predicates, `CensusSchurExclusionOffFitOn` and `CensusTorqueAliveOffFitOn`.  Three named holes remain in the ROUTE, none of them an enclosure defect: (i) the exact fit, where `generalJKernelTeacherSchurDetT` vanishes IDENTICALLY — and that is a THEOREM, `generalJKernelTeacherSchurDetT_eq_zero_at_exactFitSwap` (both unit-teacher values of the first eliminated angle equation vanish there, so the substituted kernel vector is `(0,0)`, and the cleared Schur determinant is homogeneous of degree four in it), so the un-restricted `CensusSchurExclusionOn` is FALSE over any region containing the fit (`not_censusSchurExclusionOn_of_exactFit_mem`) and row 7b's mass-carrying route is not an alternative there but the only route; (ii) the antipodal locus `theta0 - theta1 = pi`, where the ceiling `dmax < pi` is SHARP (`exists_opposite_mem_censusStripJ`, and worse for the labelled square, `exists_opposite_mem_censusSquareJ`) and five of the 55 certified non-exact-fit families sit EXACTLY, all on the equal-teacher-mass line `y = -1/2` (`replay_C_F4_0p3_m0p5.log` and `replay_C_F5_{2p4,2p6,2p8,3}_m0p5.log`, each carrying `D = 3.141592654`) — they are outside `GeneralJSeparatedPairCritical` by that predicate's own non-opposition FIELD, so no precision reaches them; (iii) the torque conditions, which the tree can only supply region-wise under both-masses-positive -- now REDUCED rather than merely named: `censusMassesJ_pos_of_mem_Ioo` supplies the positivity the tree lacked (it had only `censusMassesJ_ne_zero`), and `censusTorqueAliveOffFitOn_of_studentsOffCells` trades the torque enclosures for the geometric `CensusStudentsOffTorqueCellsOn` -- each student off the two open cells `(0,beta)`, `(pi,pi+beta)` modulo a full turn.  Stated, not discharged |
| **the dependency law that predicts which row-7a families are reachable** | measured (mean-value interval enclosure width against a dense-sampled true range, `row7a_stagewidth.py 160 all` -> `row7a_depfactor.log`; required box sides from `row7a_band.log`) | the 55 non-exact-fit families certified at the replayed witnesses, at box side `1e-4`, 160-bit mean-value form of `generalJKernelTeacherSchurDetT` with `beta` swept as an interval of the same width | the dependency factor `F` = (enclosure width) / (true range) obeys **`F ~ D^-9.17`** over four decades of student gap, with 2.3 decades of residual scatter — so the gap predicts the enclosure on average but not exactly (a family at `D = 1.4530` carries `F = 1.4e3`).  `F` in turn orders the families as the measured required box side does: Spearman **`-0.910`** on the 51 families the two tables share.  The rule that falls out is one-sided and exception-free on this data: **`F >= 1e3` REFUSES a family** — 0 of 15 is resolvable at a box side `>= 1e-4`, at any threshold from `1e3` to `1e6` — while `F < 1e3` is right 32 times in 36, the four failures all sitting in the top decade below the threshold and all missing by exactly one decade.  Scope: the determinant test only, these families only, and the sampled true range UNDER-estimates a range, so `F` is over-estimated — the direction that makes a refusal conservative rather than optimistic |
| **row-7a exclusion certificates on named gap bands** | **one valid map-only smoke; five historical SAFE artifacts invalidated** | `row7a_smoke` still replays on `[0.7500,0.7510] x [-0.4000,-0.3990]`, gap `[2.00,2.10]`, because all 26 leaves are MAP exclusions.  `row7a_safesmoke` and the four `row7a_scale*` artifacts are retained only as regression fixtures | the mass-free evaluator formerly used `sin(p)` instead of `e*sin(p)` for `|sin|` on odd half-branches.  After correction, every stored SAFE leaf in the first four invalid artifacts fails replay (142/152/229/560 mismatches), and the largest artifact fails for the same reason; none licenses a Schur exclusion.  The corrected literal evaluator, its AD derivative, and an independent canonical-`abs(sin)` gate now agree.  The selected-minimum route below supersedes this all-zero route for saddle folds |
| **selected-minimum F4 rectangle** | **certified count = 1** (`minimum_dfs_f4probe_p000.json`, four `minimum_boundary_f4probe_*` artifacts, and `minimum_witness_f4probe.json`, all independently replayed) | every `(beta,y)` in `[0.7,0.8] x [-0.45,-0.35]`, on the physical strip `censusStripJ 0.137 0.001 3.13`: the regularity tree has 94,747 leaves (83,861 map exclusions, 100 negative selectors, 10,786 positive determinants, zero undecided); the four boundary trees have 1,403 leaves; the witness `(0.75,-0.4)` has exactly one selected class, the exact fit | this earlier certificate closes one rectangle assigned the F4 label, not a global nonrectangular cell. The witness replay joins a 34-leaf map-only coincidence collar to the existing 59,233-leaf exhaustive enumeration, filters its fold copy above `D=3.13`, and evaluates Cramer `T00/det` signs directly |
| **full-zero 68-tile atlas** | externally replayed interval obligations plus the conditional Lean atlas propagation theorem; full zero count `N = 2` and selected-minimum count `N = 1` | every teacher in `[0.5,2.2] x [-0.45,-0.1]` on `censusStripJ 0.137 0.001 3.13`; 68 regularity forests: 41,371,573 leaves and 206,853,513 nodes; 272 boundary forests add 41,724 leaves and 66,040 nodes; all-forest total 41,413,297 leaves and 206,919,553 nodes, with zero undecided leaves and zero unresolved volume | the checker-to-Lean correspondence for `censusAngleMapJ`, `censusCramerSchurT00J`, and `censusCramerSchurDetTJ` is human-checked, not a Lean theorem; the concrete tile constants are external data |
| separated mass-free elimination (the `separatedDet` argument) | proved for the CENTERED kernel (`SeparatedBoundary.lean`); re-derived symbolically + numerically for this kernel | — | the noncentered instance is not itself a Lean statement; the tree's noncentered elimination is the second-order (Schur) one |


## The separated-family row, in full

The table's separated-family row is a one-line claim with three episodes
behind it.  They are kept here so the row stays scannable and nothing is lost.

**The Krawczyk Jacobian sign error (2026-08-08, commit `5e5adf1`).**  In
`census_cert.sep_res` the interval `dG0/d(th1)` was an INVALID enclosure — it
missed the true derivative range in 188 of 300 random boxes — which breaks the
Krawczyk hypothesis and voided every verdict computed with it.  Fixed, guarded
permanently by `sep_jacobian_audit.py` (a FALSIFIER: finite differences can
refute an enclosure, never certify one), and every stored point re-run item by
item (`jacobian_fix_delta.json`: 301 points, 26 changed).  The RESIDUAL itself
was always correct — it agrees with the independent float mirror to 9 digits —
so the locator, the whole coincident half, and every Lean theorem were
untouched.

**The count-matching route (2026-08-08, `count_match.py`).**  The Lean tree
holds the reduction layer (`Planar/SignedN2/FreeMassJ/SeparatedCount/`): the
separated matrix is ONE scalar `separatedNumJ` at four arguments, each student
of a family sits at a zero of a single 1-D load, and the mass-free
elimination, sign separation, lattice zeros, branch form and antipodal-slice
positivity are all axiom-clean.  The scalar's strict SIGN LAW is proved
(`SignLawJ.lean`), with the positive-mass STRADDLE (`LatticeStraddleJ.lean`)
and the mixed-sign SAME-SIDE half (`MixedSameSideJ.lean`) completing the
geometry over every teacher with `s0*s1 != 0`.  Note which sector is which:
`(s0,s1) = (sin psi, cos psi)`, `psi = (y+1)pi/2`, so `y < 0` (F2–F7) is the
POSITIVE-mass sector and `y > 0` (F1, collar) the mixed one — `count_match.py`
had these labels inverted before that pass.  What this licenses is an
ORDER-ONLY box exclusion: four endpoint comparisons, no interval evaluation of
the residual.  What it does NOT give is a proved bound on the family count;
`count_match.json` records every proposed label as `open` and none is closed by
that route.

**Two negative results, kept because each closes a route.**
(a) The arc localisation does NOT reduce the family count: the straddle /
same-side condition is IMPLIED by the load equation, so every zero of the 1-D
load already lies in the admissible arcs (measured at all 30 face witnesses)
and up to 3 bunch in ONE arc.  (b) The ambient 4-dimensional space
`span{W(.,D), W(.-beta,D), hJ, hJ(.-beta)}` is NOT a weak Chebyshev system — an
explicit element has 6 sign changes per period, verified at 24k/96k/400k/2M
grid points and at 50 and 60 decimal digits — which kills the
collocation-determinant / Dodgson route to a branch cap.  The witness is far
off the teacher-consistency variety (`c0c3 - c1c2 = -0.134` at `|c| = 1`), so
the branch cap for the student load itself is untouched; what is established is
that the coupling `(c0,c1) || (c2,c3)` is ESSENTIAL and no argument treating
the four generators as an unstructured space can prove the cap.

**Where the branch and bound stands.**  The order-only prefilter is what makes
the per-teacher enumeration tractable: it discards 15.89 of F4's 19.61 domain
area in 2 901 kills, and 2.64 of F1's 19.93.  With it, and with the
division-free centered form, the enumeration now runs to COMPLETION at
witnesses carrying all seven proposed labels (`ivcert_*.bits.json`, replayed by
`ivcert_check.py`), F1's included.  The prefilter is genuinely weaker in the
mixed sector — 39 % of the area at the F1 witness `(2.60, 0.20)` against 81 %
in the positive one — which is why F1 was last.

## What no instrument here can see

* The mass-free (row 7a) route at an exactly antipodal student pair
  `theta0 - theta1 = pi`, and at any student gap below about `0.4`.  The first
  is definitional -- non-opposition is a FIELD of
  `GeneralJSeparatedPairCritical` -- and five certified families sit there.
  The second is an enclosure limit with a measured law behind it: the
  dependency factor `F ~ D^-9.17`, and `F >= 1e3` refuses the family at any
  box side tried.  In both cases what the census actually uses there is the
  mass-CARRYING route, not this one.
* A separated family whose student gap is below `1e-4`, or one living
  inside the uncertified strips above.  (The collar families' gaps scale
  linearly with the distance to the stratum, so the `D < 1e-4` blind spot
  corresponds to `y` within about `2e-4` of the strata.)
* A second-order selector cannot decide a Schur-degenerate point by itself.
  The public two-student classification theorems give exact variational labels;
  an external interval selector still needs separate evidence when it reports
  a concrete nondegenerate sign.
* Anything about kernels, widths, or teacher counts other than the two
  models and `n = m = 2`.
