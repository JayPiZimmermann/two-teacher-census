# Cluster agent brief: close every `sorry` in `TwoTeacherCensus/`

You are working in the Lean 4 repository `two-teacher-census` (Lean
`leanprover/lean4:v4.4.0`, Mathlib pinned at `v4.4.0` through
`lake-manifest.json`; a Mathlib cache is present, `lake build` takes seconds per
changed file).  Read `PROOF_PLAN.md` first: it lists, theorem by theorem, the
lemma chain, the manuscript Claim/Step each lemma realises, the expected
difficulty, and the long-tree declarations whose *mathematics* (not code) you
may consult in
`/home/jzimmermann/code/Skip-Connections-Avoid-Spurious-Local-Minima/LeanFormalization`
if that path is available; the paper repository
`Removing-spurious-minima-for-planar-features-by-skip-connections` holds the
Gaussian pair-moment calculus (`Preliminaries.lean`) that
`L_eq_gaussianLoss` needs.

## Goal

Make every `sorry` in `TwoTeacherCensus/Definitions.lean`,
`TwoTeacherCensus/Classification.lean` and `TwoTeacherCensus/CensusMap.lean`
disappear, so that `lake build` reports 0 errors and 0 `declaration uses
'sorry'` warnings and `cluster/audit_axioms.sh` reports only `propext`,
`Classical.choice`, `Quot.sound`.

Exception: `plain_census_map_conjecture` in `CensusMap.lean` is a conjecture
whose docstring says so; leave its `sorry` and do not spend effort on it.
Everything else is a theorem of the manuscript and is expected to be provable.

## Rules

1. **Statements are fixed.**  Do not change the statement of any theorem,
   definition or `inductive` unless it is mathematically wrong.  If you are
   convinced a statement is wrong (a missing hypothesis, a wrong constant, a
   wrong sign), fix it minimally and record the change and its justification in
   `PROOF_PLAN.md` under "Statement corrections", with the manuscript equation
   or long-tree declaration that decides it.
2. **No new axioms.**  No `axiom`, no `native_decide`, no `Lean.ofReduceBool`,
   no `sorry` in anything you call finished.  `decide` on finite types is fine.
3. **No imports from the old trees.**  Only `import Mathlib` and modules of this
   repository.  You may *port* proofs from the paper repository
   (`Preliminaries.lean` etc.) by copying them into this tree with the
   `TwoTeacherCensus` namespace; cite the source in a comment.
4. **Prefer shorter proofs.**  The point of this repository is a development of
   a few files and a few thousand lines.  Factor shared calculations into
   lemmas with docstrings; do not duplicate 100-line proofs for the two models
   when a `Model`-indexed lemma works.
5. **Run `lake build` after every change** and fix errors before moving on.  A
   file with errors is worse than a file with `sorry`.  Never leave the build
   broken at the end of your pass.
6. **Keep docstrings naming the paper Claim/Step** (`clm:census-*`,
   `stp:census-*`, `prop-two-student-system`, `prop-collision-block`,
   `lem-flat-sufficiency`, `thm-map-centered-complete`, ...).  When you add a
   lemma, say which Step it serves.
7. **Order of work.**  Follow the dependency order in `PROOF_PLAN.md`:
   first the bridges everything rests on (`isCritical_iff_fourRows`,
   `differentiableAt_lossPair`, `Δ_pos_of_not_modEq`, `Δ_eq_zero_iff`,
   `L_eq_gaussianLoss`), then the collision block
   (`collided_isLocalMin_iff_selector`, `collided_not_isGlobalMin`), the dead
   students (`dead_not_isLocalMin`, three-translate rigidity), the per-stratum
   classification lemmas, the per-family type lemmas, and last the separated
   second-variation results (`centered_separatedBeam_saddle`,
   `plain_separatedSmooth_type`, `centered_separatedBeam_existsUnique`) and
   the census map.
8. **One pass = a few lemmas.**  Pick the lowest-level unproved lemmas, prove
   them, build, and stop.  The driver script re-invokes you with the new build
   state; do not try to do everything in one invocation.
9. **Do not edit** `cluster/`, `web/`, `certificates/`, `deploy/`, `README.md`
   or `lakefile.lean` (except `roots` if you add a module).  Do not commit; the
   driver commits.
10. If you add a module, add it to the `roots` list in `lakefile.lean` and
    import it from `TwoTeacherCensus.lean`, and keep the total under five
    modules.

## Stopping condition

You are finished for this pass when `lake build` succeeds with strictly fewer
`sorry` than before and no errors.  The whole job is finished when
`lake build` has 0 errors, 0 `sorry` warnings outside
`plain_census_map_conjecture`, and `cluster/audit_axioms.sh` is clean.  Then
write a short "Status" paragraph at the top of `PROOF_PLAN.md` and stop.
