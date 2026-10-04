# Two-teacher census

The complete classification of the critical points of a planar two-student
network trained against a two-atom teacher, for the centered feature
`C(z) = |z|/2` and the plain ReLU `R(z) = max(0, z)`, under the Gaussian
population loss — every critical point, each with its proved variational
type (global minimum, spurious local minimum, topological saddle) — together
with the **census map**: for which teachers the list of minima is the same.
The mathematics is Section 5 and Appendix C of the paper
[arXiv:2610.01728](https://arxiv.org/abs/2610.01728); the original (long)
formalization is the `LeanFormalization` tree of the research repository
`Skip-Connections-Avoid-Spurious-Local-Minima`.  This repository is the short,
self-contained re-proof and the public explorer.

Interactive explorer: <https://two-teacher-census.fly.dev> (placeholder until
the first deploy).

## Lean development (`TwoTeacherCensus/`)

Status: **statements fixed, proofs in progress on the cluster.**

| module | content |
|---|---|
| `Definitions.lean` | the two kernels `Φ_C, Φ_R` and couplings `H_C, H_R` (verbatim from the paper repository), the planar losses `L_C`, `L_R` in kernel form, criticality (`fderiv = 0`), local/global minima, topological saddles, the two-atom teacher data `P_q, A_q, W, τ_q`, Cramér masses, canonical gaps, and the census; all 51 lemmas proved (kernel calculus, parity, periods, branch formulas) |
| `Classification.lean` | the two census theorems with every table row explicit: `centered_two_student_classification : IsCriticalL_C β s c θ ↔ CenteredFamilies β s c θ` (ten rows), `centered_two_student_type`, and the plain-ReLU analogues `plain_two_student_classification` (fourteen rows), `plain_two_student_type`; the lemma chain of Appendix C0 (reduction isometries, four rows, `L_R = L_C + ⅛‖u−v‖²`, collision block, dead-student revival, Schur type) with `sorry` where the proof is pending |
| `CensusMap.lean` | the map coordinates, the three mass-free determinants and the sign bridge (proved), the centered census-map theorem `centered_census_map` (nine curves, ten faces, four values, five regions), and the plain-ReLU map as a documented conjecture |

`PROOF_PLAN.md` gives, theorem by theorem, the lemma chain, the long-tree
results each lemma replaces, the estimated difficulty and the open risks.
The proofs are being completed by an automated agent loop
(`cluster/perfectionize.sh` driving `cluster/AGENT.md`); `cluster/audit_axioms.sh`
checks that nothing beyond `propext`, `Classical.choice`, `Quot.sound` is
used, and `cluster/build_certificate.sh` writes `BUILD_CERTIFICATE.txt`.

## Certificates (`certificates/`)

Interval-arithmetic and symbolic certificates for the census map, copied from
the research repository: `census_map_centered/` (root counts, boundary
curves, per-face witnesses), `census_map_noncentered/` (coincident counts,
Krawczyk-certified families, continuation and atlas certificates).
`COMPLETENESS.md` states, claim by claim, what each instrument establishes and
what it cannot see; files above 3 MB are omitted and listed with hashes in
`LARGE_ARTIFACTS.md`.

## Explorer (`web/`)

A static page: `index.html` with the classification tables, `explorer.js`
(the landscape explorer widget), and the precomputed `census-map.js`, which
`web/precompute/precompute_census_map.js` generates from the classifier in
`explorer.js` (`node precompute_census_map.js --depth 8`; `--check` verifies
the shipped asset).  `web/PROVENANCE.md` records the source commit of every
file.  `deploy/` and `fly.toml` host it on fly.io.

## How to build

```bash
# Lean (needs the Mathlib cache; a few seconds per module afterwards)
lake exe cache get
lake build
cluster/audit_axioms.sh          # axiom audit of every theorem

# explorer
python3 -m http.server --directory web 8765
```

Toolchain `leanprover/lean4:v4.4.0`, Mathlib pinned in `lake-manifest.json`
(the same pins as the paper repository).
