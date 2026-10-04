# Proof plan: the two-teacher census in five files

Status (2026-10-03, laptop draft): statements fixed; 168 theorems in three
modules; 98 axiom-clean (checked with `#print axioms`), 70 depend on `sorryAx`
(48 carry a literal `sorry`, 22 are proved modulo sorried dependencies).
Elaboration of every module checked with `lean` against the pinned Mathlib
(`v4.4.0` toolchain and manifest copied from the paper repository); no
`lake build` has been run yet on this tree.

This document is the architecture the cluster agent (`cluster/AGENT.md`)
follows.  For each theorem: the lemma chain, the long-tree results it
replaces, the estimated difficulty, and what remains to do.  Difficulty scale:
**S** (an afternoon: algebra, case splits, direct Mathlib API), **M** (a
day: real analysis with explicit inequalities, second-variation arguments),
**L** (several days: the separated stratum, Gaussian pair moments, root
counts), **XL** (open or certificate-only).

The long tree being replaced lives at
`Skip-Connections-Avoid-Spurious-Local-Minima/LeanFormalization` (843 modules;
public surface `Main/TwoStudentClassification.lean`, paper route
`Main/Paper/Appendix/C0_TwoTeacherCensus.lean`).  Nothing is imported from it;
its declarations are cited as the mathematical reference only.

## 0. Design decisions that shorten the proof

* **One `Model`-indexed development.**  `Φ q`, `Hq q`, `κ q`, `period q`
  carry both kernels; every lemma that the long tree proves twice
  (`...Centered` / `...J`) is stated once.  Only the strata differ: the
  centered reduction has a projective gap `[0, π)`, the plain one an oriented
  gap `[0, π]` with the extra antipodal stratum.
* **Frame-free families.**  The rows of the tables are stated in the original
  frame with congruences `≡ [PMOD period q]` (the bisector `(β₀+β₁)/2`, the
  interlacing predicate via `dirRep`), so no `centeredTranslate`/`dirRep`
  transport layer between "canonical" and "original" labels is needed; the
  translation isometry is proved once (`isCritical_translate_iff`,
  `isLocalMinL_translate_iff`, `isGlobalMin_translate_iff`) for the proofs
  that want to move the first teacher to `0`.
* **Split-free collision rows.**  The collision line is stated with free
  split in every stratum (the long tree's mixed stratum lists collided
  dead students separately); the selector handles `c₀c₁ = 0` correctly, and
  the dead rows overlap harmlessly.
* **The "real corner" row is a theorem**, `plain_realCorner_forces_exactFit`,
  not a family: its content is that it coincides with the exact fit.
* **The census has multiplicity** (`CensusValue`), so `fit + coin⁻` and
  `fit + 2 coin⁻` are different values; `Census` as a set is kept for the
  `(kind, type)` legend.
* **`L` is the finite kernel energy** (the manuscript's `eq-loss-quadratic`);
  the literal Gaussian integral `GaussianLoss` appears only in
  `L_eq_gaussianLoss`, which is the single analytic input (`L ≥ 0`).  All
  zero-loss (global) rows reduce to algebra (`L_eq_zero_of_exactFit`,
  `L_eq_zero_of_zeroMeasure`, proved) plus this one lemma.

## 1. `centered_two_student_classification` (`thm:centered-two-teacher-census`)

`IsCriticalL_C β s c θ ↔ ∃ f, CenteredMember f β s c θ` — **assembly proved**
from three stratum lemmas by `centered_strata_exhaustive` (proved).

| lemma | realises | replaces (long tree) | difficulty | status |
|---|---|---|---|---|
| `isCritical_iff_fourRows` | `prop-two-student-system` | `isCriticalL_C_two_student_iff_four_equations`, `freeMassGeneral_system_iff` | M | sorry |
| `differentiableAt_lossPair` | (used by the above) | `kernelEnergy_differentiable` (paper repo) | S | sorry |
| `Δ_pos_of_not_modEq`, `Δ_eq_zero_iff` | `prop-kernel-properties`(6), `eq-classification-radial-determinant` | `massBlock_det_pos`, `phiCos_lt_pi_div_two_of_...` | M | sorry |
| `radialRows_iff_cramer` | `eq-separated-cramer` | `centeredMasses_eq_cramer...` | S | **proved** |
| `angularRows_iff_torqueMasses` | `eq-classification-angular-derivative-masses` | — | S | **proved** |
| `zeroField_critical_iff` | zero-field rows | `centered_zeroTeacher_critical_iff_exact_global` | S | **proved** (modulo rows, Δ) |
| `centered_oneLine_critical_iff` | `stp:census-centered-one-line` | `centered_oneTeacher_critical_iff_exact_global_or_perpendicular_saddle` | M | sorry |
| `centered_genuine_critical_iff` | `stp:census-centered-zero-mass`, `-quarter`, `-beam`, `-mixed` | `centered_positiveTeacher_critical_iff_complete_classification`, `centered_negativeTeacher_...`, `centered_mixedTeacher_critical_iff_geometric_classification` | L | sorry |

Plan for `isCritical_iff_fourRows` (the bridge everything rests on):
`lossPair` is a finite sum of `c i * c j * Φ q (θ i − θ j)` terms, each
differentiable on the product (`hasDerivAt_Φ`, `DifferentiableAt.comp` with
the coordinate projections `(differentiable_apply i).comp differentiable_fst`).
Then `fderiv = 0 ↔ ∀ basis vector v, fderiv p v = 0`, and
`fderiv p (Pi.single i 1, 0)` is the derivative of the one-variable mass
curve (`HasFDerivAt.comp_hasDerivAt` along `r ↦ (update c i (c i + r), θ)`),
which `L_mass0_curve` (proved) already identifies with `R_0/(2π)`; the
angle directions use `hasDerivAt_Φ`.  About 120 lines.

Plan for `Δ_pos_of_not_modEq`: the window `1 ≤ Φ_C ≤ π/2`, `0 ≤ Φ_R ≤ π`
with equality cases.  Write `Φ_C(t) = g(|cos t|)` with
`g(u) = u arcsin u + √(1−u²)`, `g' = arcsin ≥ 0` on `[0,1]`
(`hasDerivAt_orthantPhi`, proved) so `g ≤ g(1) = π/2` with equality iff
`|cos t| = 1`; for `Φ_R` add `(π/2) cos t` and use `Φ_R ≥ 0` from the branch
formula (`phiCosJ_branch`, proved) and `tan u ≥ u` (`Real.lt_tan`).

Plan for `centered_oneLine_critical_iff`: after `isCritical_translate_iff`
put the line at `0`; `P_C = σ Φ_C`, `A_C = σ H_C`, `H_C` vanishes exactly on
`(π/2)ℤ` (`couplingH_branch`, proved, plus periodicity).  Case split on dead
students via `fourRows_collided_iff` (proved) and `radialRows_iff_cramer`.
The exclusion of a separated two-live point is the manuscript's trick: place
the zero-mass second teacher on student 0 and use the mixed-stratum endpoint
argument — or directly: `A_C(θ₀) = A_C(θ₁) = 0` puts both students on
`{0, π/2}`, and the Cramér masses at `(0, π/2)` give `c = (σ, 0)`, dead.

Plan for `centered_genuine_critical_iff` (the long one): by student precedence
— both dead (`BothDead`, from the radial rows with `c = 0`), one dead
(`OneDead`: `fourRows` with `c i = 0` is literally the predicate), collided
(`fourRows_collided_iff`, proved), separated live.  The separated live case
is the only hard one:
* `H_C(D) = 0` with `D ≢ 0` means `D ≡ π/2`: the quarter case.  Prove
  `lem-quarter-rigidity` (`genTorque_bisector_eq`, `bisector_fold_rigid`,
  `chord_sin_gt` in the long tree): a quarter-separated pair of torque roots
  off the teacher lines forces `s₀ = s₁` and the bisector pair, or sits on
  the teacher lines with `β̂ = π/2` (the exact fits).  Difficulty M (one
  monotonicity, `Ξ(x) = (π/2 − x) tan x / x` strictly decreasing via
  `sin y > y(π−y)/π`).
* `H_C(D) ≠ 0`: `TorqueMasses` and `CramerMasses` both hold; for same-sign
  teachers the balance-scalar sign law (`clm:census-beam-signs`,
  `separatedUnitBalanceNum_pos/neg`) forces interlacing and `s₀ ≠ s₁`
  unless at the exact fit or the quarter pair (`generalCenteredSeparated_nonExact_iff_beamSaddle`); for
  mixed teachers the same sign law forces the exact fit
  (`stp:census-centered-mixed`, `prop-map-mixed-rigidity`).  Difficulty L:
  this is the content of `Planar/SignedMass/SignedTeacher/Classification/
  {SameSign*,Mixed*,Beam*}` (~40 modules); the manuscript's balance-scalar
  route (`eq-map-centered-balance-scalar`, `eq-map-centered-scalar-sign`) is
  the short version and is what to formalise.

## 2. `centered_two_student_type` (type clause)

**Assembly proved** from ten per-family lemmas.  `oneDead`, `bothDead` are
proved from `dead_not_isLocalMin` + `isTopologicalSaddle_iff_not_isLocalMin`
(proved); `zeroMeasure`, `exactFit` are proved from `L_eq_zero_of_*` (proved)
+ `L_nonneg` (needs `L_eq_gaussianLoss`).

| lemma | realises | replaces | difficulty | status |
|---|---|---|---|---|
| `L_eq_gaussianLoss` → `L_nonneg` | `eq-loss-quadratic`, `lem:pair-moments` | paper repo `eq_loss_quadratic` (~1300 lines of `Preliminaries.lean`: rotation invariance, polar pair moment `E|⟨u,x⟩||⟨z,x⟩| = (2/π)(ρ arcsin ρ + √(1−ρ²))`) | L (port) | sorry |
| `centered_oneLineExact_global` | one-line exact rows | `centered_oneTeacher_...exact_global` | S (algebra + `L_nonneg`) | sorry |
| `centered_oneLinePerp_saddle` | perpendicular row | `...perpendicular_saddle` | M: rotation formula `L_collided` (proved) gives `2πL(t+r) = κμ²/2 − μ σΦ_C(π/2 + r) + E/2`, second derivative `−(2/π)σ² < 0` via `hasDerivAt_couplingH` | sorry |
| `centered_collision_type` | `stp:census-centered-collision` | `SignedCenteredLineCritical.exact_variational_label`, `MixedLineMinCondition` | from `collided_isLocalMin_iff_selector` + `collided_not_isGlobalMin` | sorry |
| `collided_isLocalMin_iff_selector` | `prop-collision-block`(4), `clm:census-centered-selector`, `lem-flat-sufficiency`, `stp:census-flat-*` | `coincidenceAngularBlock_det_trace_pos_iff`, `centered_flatCollision_localMin_iff_paper`, `collisionQuadratic_positiveDefinite_iff_selector` | L: exact decomposition `2πΔL = (κ/2)S² − âb̂δ(x−y) + âΔP(x) + b̂ΔP(y)` (algebra, like `L_revive_dead`), second-order Taylor of `P_q` with Lipschitz `P''` (`hasDerivAt_Hq`), definiteness ⇔ selector (`collisionAngularBlock_det_trace`, proved), flat case via the quartic profile `P(t+r) = P(cos r + (r/2) sin r)` (ODE uniqueness on a cell) and the floor `(κ/2)ϖ² + (|ab|κ/16)(x−y)² + (μP/192)y⁴` | sorry |
| `collided_not_isGlobalMin` | `prop-circle-l2` | `centered_collision_loss_pos` | M: `L_collided` (proved) gives `2πL = E/2 − P(t)²/(2κ)` at the pin; need `E κ > P(t)²` for two genuine lines, a Cauchy–Schwarz in the kernel (or: `L = 0` forces the exact fit, by `L_eq_gaussianLoss` and the three-translate rigidity) | sorry |
| `centered_quarterPair_saddle`, `centered_orthogonalDiagonal_saddle` | `stp:census-centered-quarter` | `quarter_strict_saddle`, `GenQuarterBranch.exactGlobal_or_nonExact_not_isLocalMin` | M: explicit negative direction of the 4×4 Hessian with mass block `[[π/2,1],[1,π/2]]` | sorry |
| `centered_separatedBeam_saddle` | `stp:census-centered-beam` | `CenteredPositiveBeamSaddle`, `Beam.*` (~25 modules) | **L — the hard part**: second variation at the beam root; the long tree certifies descent along an explicit direction using `clm:census-beam-signs` (`beam_signs`, stated) | sorry |
| `centered_separatedBeam_existsUnique` | structure clause | `centeredBeamRoot_existsUnique`, `centeredBeamScale_strictAntiOn` | L: phase defect `F_{β̂,a}(b)` strictly decreasing with endpoint signs, scale `S_{β̂}` strictly decreasing from `+∞` to `0` | sorry |
| `beam_signs` | `clm:census-beam-signs` | `Beam.beamSigns_paper` | M: two trigonometric inequalities on the beam triangle | sorry |

## 3. `plain_two_student_classification` (`thm:plain-two-teacher-census`)

**Assembly proved** from four stratum lemmas by `plain_strata_exhaustive`
(proved).

| lemma | realises | replaces | difficulty | status |
|---|---|---|---|---|
| `plain_zeroField_critical_iff` | zero-field rows | `coincidentTeacherJ_critical_iff_labeledClassification` (zero part) | S | **proved** (modulo rows, Δ) |
| `plain_oneRay_critical_iff` | `stp:census-plain-reduction` one-ray rows | `atMostOneTeacherJ_critical_iff_labeledClassification`, `secondOnlyTeacherJ_...` | M: `H_R` vanishes exactly on `πℤ`, `Φ_R(π) = 0` (`phiCosJ_pi`, proved), antipodal cancellation `c₀ + c₁ = 0` at `a + π` | sorry |
| `plain_antipodal_critical_iff` | `stp:census-plain-antipodal-collision`, `-separated` | `antipodalTeacher_critical_iff_canonicalLabeledClassification`, `AntipodalPaperStepSurfaceJ` | L: torque factorisation `A_R(t) = [πs₀ − St] sin t` on `(0,π)` (from `couplingHJ_branch`, proved), roots `0, π, t±`; separated families via the midpoint equations | sorry |
| `plain_genuine_critical_iff` | `stp:census-plain-gap-collision`, `-separated`, `clm:census-plain-load-sign`, `clm:census-plain-ghost-diagonal` | `generalJ_genuineTeacher_critical_iff_labeledClassification` (~60 modules under `SignedN2/FreeMassJ`) | L: the exhaustive part is the student-precedence case split (as in §1); the separated live case is *by definition* `CramerMasses ∧ TorqueMasses` plus the ghost/smooth/corner trichotomy, and the corner case needs `plain_realCorner_forces_exactFit` | sorry |
| `plain_realCorner_forces_exactFit` | real-corner row | `generalJ_realCorner_exactFit` | M: a student on a teacher ray makes the one-sided third derivatives of `P_R` differ by the teacher mass (`eq-flat-boundary-derivatives`), incompatible with the smooth torque row unless the other student completes the fit | sorry |
| `plain_load_sign` | `clm:census-plain-load-sign` | `SeparatedCount/SignLawJ.lean` | M | sorry |

## 4. `plain_two_student_type`

**Assembly proved**; `oneDead`, `bothDead`, `zeroMeasure`, `exactFit` as in §2.

| lemma | realises | replaces | difficulty | status |
|---|---|---|---|---|
| `plain_oneRayExact_global` | one-ray exact | — | S | sorry |
| `plain_oneRayAntipodalCancel_saddle` | `stp:census-plain-reduction` | `PlainAntipodalCancelPaper` | S/M: revive at `x` near `a+π`: `L_revive_dead` (proved) with `F(x) = −σΦ_R(x−a) < 0` | sorry |
| `plain_antipodalInteriorCollision_saddle` | `stp:census-plain-antipodal-collision` | `antipodalInteriorCollision_saddle` | M: selector with `τ_R` at `t±` | sorry |
| `plain_antipodalEndpointCollision_type` | endpoint cusp selector `N > 0` | `antipodalEndpointNullCubicJ`, `AntipodalEndpointNullCubicClassificationJ` | L: the second variation is degenerate along one direction; the cubic `|z|³/3` corners of the students versus the teacher's corner decide (`eq-classification-antipodal-expansion`) | sorry |
| `plain_oppositePerp_saddle` | opposite perpendicular | — | M | sorry |
| `plain_separatedMidpointPair_saddle` | `stp:census-plain-antipodal-separated` | `AntipodalSeparatedSymmetricReducedJ.isSaddle` | M/L: Schur bracket `generalJCanonicalSchurBracket_eq_load_deficit` | sorry |
| `plain_collisionRay_type` | `stp:census-plain-gap-collision` | `noncentered_flatCollision_localMin_iff_paper` | from `collided_isLocalMin_iff_selector` (`q = .plainRelu`) | sorry |
| `plain_oppositeBisector_saddle` | opposite bisector | `GeneralJOppositeBisector.saddle` | M | sorry |
| `plain_separatedGhost_saddle` | `clm:census-plain-ghost-diagonal` | `SeparatedGhostDiagonalJ` | M: `𝒜_ii < 0` at the incident student, so the raw angular block has a negative diagonal entry; descent along that coordinate | sorry |
| `plain_separatedSmooth_type` | `stp:census-plain-gap-separated`, `thm-separated-index-law` | `GeneralJSeparatedPairCritical.*`, `generalJSchurT00_eq_neg_c0_jac00`, … | **L — the hard part**: 4×4 second variation, Schur complement `T` (`schurT00/11/01`, defined), `T ≻ 0 ⇒` strict local min (Taylor with the `C²` kernel), negative direction `⇒` descent; the rank-one case is stated only as `Spurious ∨ Saddle`, which follows from non-globality (`collided_not_isGlobalMin`-style) and `variationalType_of_isCritical` | sorry |

## 5. Shared analytic lemmas

| lemma | realises | difficulty | status |
|---|---|---|---|
| `dead_not_isLocalMin` | `stp:census-centered-revival`, `clm:census-signed-revival` | M given rigidity: `L_revive_dead` (proved) + a direction `x` near `θ_i` with `F(x) ≠ 0` (else `F ≡ 0` on an arc, contradicting three-translate rigidity) | sorry |
| `centered_three_translates_rigidity`, `plain_three_translates_rigidity` | `clm:census-centered-three-translates`, `clm:census-plain-three-translates` | M: on an arc avoiding the kinks the three translates are trigonometric polynomials `α cos + β sin + γ t cos + δ t sin`; the beam minor `(b−a) sin a sin b` is nonzero for distinct lines (`threeKernel_zeroArc_weights_zero`) | sorry |
| `W_ne_zero_of_torqueRoot_*` | `clm:census-collision-weight` | M | sorry |
| `τ_ne_zero_of_torqueRoot_mixed` | `clm:census-centered-mixed-curvature` | M | sorry |
| `isCritical_both_iff_matched` | `one-harmonic-apart` | S/M from `L_R_eq_L_C_add_firstMoment` (proved) and the row corrections `F_R = F_C + (π/2)⟨w, e⟩` | sorry |

## 6. `CensusMap.lean`

| lemma | realises | difficulty | status |
|---|---|---|---|
| `wTau_add_wPot`, `massLine_of_torqueRoot`, `signBridge`, `boundary_elimination` | `eq-map-wronskian-relation`, `eq-map-mass-line`, `eq-map-sign-bridge`, `prop-map-boundary-elimination` (both models) | S | **proved** |
| `existsUnique_betaStarC2`, `existsUnique_uStarR` | `eq-map-centered-betastar`, `eq-map-plain-relu-betastar` | S/M (monotone function on an interval) | sorry |
| `centered_root_counts` | `thm-map-centered-counts` | L: the reduced forms `Υ ± Ω` in `(b, w)` coordinates, five inequalities C1–C5 (certified numerically in the long tree; a Lean proof needs the explicit trigonometric estimates) | sorry |
| `centered_curves_structure` | lens is a Jordan curve over `[β*₁, β*₂]`, potential curve two branches | L | sorry |
| `centered_census_map` | `thm-map-centered-complete` | L: from `centered_two_student_classification`/`_type` (which rows are minima), `signBridge` (selector in terms of determinants), `centered_root_counts` (how many roots), alternation of maxima/minima of `P_C`, and the sign of `P_C` at the roots in the mixed sector | sorry |
| `centered_census_locally_constant` | `thm-param-zero-count` | M (implicit function theorem for simple roots) | sorry |
| `plain_root_counts`, `plain_colliding_census_of_signChart` | `thm-map-plain-relu-counts`, `eq-map-selector-signs` | L / S | sorry |
| `plain_census_map_conjecture` | plain map | **XL — conjecture**, stays `sorry` (allowlisted) | sorry |

## 7. Open mathematical risks (honest list)

1. **Separated stratum second variation** (`centered_separatedBeam_saddle`,
   `plain_separatedSmooth_type`, `centered_separatedBeam_existsUnique`).  In
   the long tree this is ~60 modules of explicit Hessian algebra and two
   monotonicity arguments.  The manuscript's route (balance scalar, Schur
   complement as the balance map's Jacobian, `thm-separated-index-law`) is
   shorter but still the bulk of the remaining work.  Nothing here is known
   to be false, but the Lean cost is the largest.
2. **The rank-one Schur boundary** is typed only as `Spurious ∨ Saddle`.  The
   manuscript gives no closed-form selector either (it describes the
   reduced one-variable Taylor test).  If a sharper statement is wanted it
   has to be designed first.
3. **Transcription risk in the antipodal midpoint pair** (`AntipodalMidpointPair`,
   from `eq-antipodal-midpoint-root`): the two scalar equations were copied
   from the manuscript table; check against
   `AntipodalSeparatedSymmetric*J` before investing in its type proof.  Same
   remark for the endpoint null cubic `nullCubic` (`eq-antipodal-endpoint-cusp-selector`).
4. **Flat collision at a kink.**  `CollisionSelector` requires `OffKink` in
   the flat branch; the manuscript says collided points on a teacher line
   with `τ = 0` are saddles (same-sign) and never flat (mixed,
   `clm:census-centered-mixed-curvature`).  If a flat kink collision turned
   out to be a minimum in some signed case the selector needs a third
   branch; the long tree's `FlatCollisionPaper*` and `MixedLineMinCondition`
   say it does not.
5. **`L_eq_gaussianLoss`** is a long but mechanical port; without it every
   "global" label rests on `sorryAx`.  Alternative: prove positive
   definiteness of the kernel directly (Fourier coefficients of `Φ_C` are
   nonnegative) — not shorter.
6. **Root counts of the map** (`centered_root_counts`) were closed in the
   long tree only by interval certificates plus exact reduced-form
   inequalities; a pure Lean proof is new work (L).  The plain map is a
   conjecture and is excluded from the completion target.

## 8. What the cluster agent must do, in order

1. `differentiableAt_lossPair`, `isCritical_iff_fourRows` (unlocks every
   classification lemma; `L_mass0_curve` and `hasDerivAt_Φ` are ready).
2. `Δ_pos_of_not_modEq`, `Δ_eq_zero_iff` (kernel window).
3. `L_eq_gaussianLoss` (port from the paper repository) — can run in parallel
   with 4–6 since only the "global" labels depend on it.
4. Collision block: exact decomposition lemma, `collided_isLocalMin_iff_selector`,
   `collided_not_isGlobalMin`, then `centered_collision_type`,
   `plain_collisionRay_type`, `centered_oneLinePerp_saddle`.
5. Dead students: three-translate rigidity, `dead_not_isLocalMin`.
6. Stratum lemmas: `centered_oneLine_critical_iff`, `plain_oneRay_critical_iff`,
   then the genuine/antipodal strata by student precedence.
7. Separated second variation (§7.1), quarter rigidity, beam structure.
8. Census map (§6), leaving `plain_census_map_conjecture`.

## Statement corrections

(none yet — the agent records any change of a fixed statement here, with
the manuscript equation or long-tree declaration that decides it)
