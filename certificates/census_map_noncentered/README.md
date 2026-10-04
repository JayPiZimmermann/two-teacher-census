# Evidence and replay guide for the plain-ReLU (noncentered) census map

**Status date: 2026-08-24.** The authoritative claim ledger is
`CERTIFICATE.md`; this file is its reproduction and engineering guide.

The companion of `../census_map_centered/`, for the noncentered kernel.  Same
question -- **are the pieces drawn in the map all of the pieces there are?** --
with a mixed answer. Exact point enumerations and one 68-tile regional atlas
are closed; the interpolated global wall arrangement is not. This model's
census has a separated half, and the shipped classifier's separated locator
was wrong before it was replaced by the enumeration developed here.

## The headline

**All seven proposed labels carry replayed witnesses.**  F1, F2, F3, F4, F5,
F6 and F7 each have at least one teacher whose completeness certificate has been REPLAYED
VALID AND COMPLETE — every leaf re-derived from the spec, zero failed, zero
undecided, class disjointness certified.  The search is not the claim; the
replay is. These are exact statements at named teachers, not proof that F1--F7
are the cells of a complete arrangement or that the count is constant on each
label.

**F1 was the last, and it closed on 2026-08-10.**  Its enumeration at the
representative teacher `(0.75, 0.40)` ran the `19.9277` domain to completion in
509 533 steps and 254 893 leaves with ZERO undecided boxes and zero undecided
area, and the replay verified all 254 893 from the packed bitstring: 6
enclosures, 1 certified merge, 5 disjointness-certified classes — **four
separated families and the exact fit**, which is the count its census
predicted.  The earlier F1 attempt had exhausted a 20 000 s budget at 4.49
million steps with a live area of `0.17591` still open; what closed it is the
division-free centered enclosure form together with the hybrid search order,
and the run that closed it took 2 370 s.

**Separately, all seven labels now also carry CONTINUATION certificates** —
nine rectangular grids, 2 589 grid teachers (2 572 distinct), 30 278 certified
steps, zero undecided leaves,
families certified pairwise disjoint at every teacher.  That is a different
pointwise claim from the enumeration and neither implies a global-cell
statement; the section on the continuation says exactly which grids and
tracked families are covered.

**The selected-minimum route now closes its first two-dimensional region.**
The F4 rectangle `[0.7,0.8] x [-0.45,-0.35]` has a replayed full-period
regularity tree, four replayed physical-boundary trees, and a replayed witness
count of exactly one selected root (the exact fit).  This proves the minimum
count is one throughout that rectangle even though an additional separated
saddle is present.  It is deliberately narrower than a claim about all F4;
the certificate accounting and reproduction commands are in
`MINIMUM_DFS.md`.

**The stricter full-zero route closes a larger 68-tile rectangle.**  On
`beta in [0.5,2.2]`, `y in [-0.45,-0.1]`, and the strip
`censusStripJ 0.137 0.001 3.13`, the external replay has 41,371,573 leaves
and 206,853,513 nodes in the 68 merged regularity forests. The 272
physical-boundary forests add 41,724 leaves and 66,040 nodes, for all-forest
totals of 41,413,297 leaves and 206,919,553 nodes, with zero undecided leaves
and zero unresolved volume.
The full zero count is two and the selected-minimum count is one. The Lean
atlas theorem is conditional and generic; the concrete tiles are external
data, and the checker-to-Lean expression correspondence is human-checked.
`MINIMUM_DFS.md` states that trust boundary in full.

| noncentered object | current status |
|---|---|
| F1--F7 and the two collar refinements | proposed global labels; witness evidence, not a certified arrangement |
| representative pointwise censuses and collar row | five exact values at the seven replayed F1--F7 representatives; the collar supplies an additional pattern whose listed family is certified to exist, not an exhaustive census |
| proposed Euler bookkeeping | `V=11, E=18, F=7`, `chi=0` on the cylinder; conditional on the interpolated arcs and absence of extra walls |
| proposed full-map doubling | `V=14, E=28, F=14`; the collar refinement to 18 is likewise conditional |
| exact teacher enumerations | twenty-two replays, 2,368,457 leaves, zero failed or undecided leaves |
| closed regional atlas | 68 tiles on the rectangle above; full count 2, selected count 1, with the stated external/Lean trust boundary |
| shipped map (commit `2201dd0`, depth 8) | 49 rendered components, 19 sub-threshold; a finite-resolution measurement |

Against the pre-fix classifier (commit `f430fab`) the comparison measured a
one-sided disagreement at 52 of 146 band points whose returned locator lists
passed their per-family interval checks, and 1868 rendered fragments; the fix
(`2201dd0`) adopted the certificate's locator as the shipped one. The current
comparison has zero disagreement wherever the locator returns an accepted
list. This is nonexhaustive returned-list agreement outside the separately
enumerated teachers, not a certified band census. `CERTIFICATE.md` section 8
records both states.

## What is structurally different from the centered model

**There is no potential curve at all.**  `Wpot` has no interior roots for
any `beta` — its only zeros are the lattice roots `t = pi` and
`t = beta + pi` — so the centered map's second boundary curve has no
analogue.  This is now a theorem-level statement of the certificate
(`certify_analytic.py`), not a box covering.

**The kink lattice has four points**, `{0, beta, pi, beta+pi}`, not two.
Closed branch forms give an exact global reduction: with `p, q` the branch
coordinates of `t` and `t - beta` and `d = q - p`,

```
   Wtau = |pq| (sin d + d sinc(p) sinc(q))
   Wpot = |pq| (−sin d + d sinc(p) sinc(q))
   Wwgt = |pq| d sinc(p) sinc(q)
```

(sympy-exact on all four sign quadrants), and `Wtau + Wpot = 2 Wwgt` is
re-derived rather than assumed (interval residual ≤ 7.5e-58).

**The coincident root counts are closed analytically, for every `beta`.**
Exact piecewise inequalities on the sinc-reduced forms — positivity of the
sinc product, `sinc(a) sinc(b) > sinc(a+b)`, and a certified unimodality
inequality on the lens piece — give, with no box covering and hence no
collars:

```
   Wtau: 2 roots outside [beta*_1, beta*_2], 3 at the endpoints (double
         root), 4 inside;   Wpot: always 2
   beta*_1 = 2u*, tan u* = pi − u*:  u* = 1.11283481547935901488567226
   beta*_1 = 2.22566963095871802977134452
   beta*_2 = 2pi − beta*_1 = 4.05751567622086844715394225
```

the exact analogue of the centered `tan u* = pi/2 − u*`.  The earlier
adaptive box coverings (309 + 305 boxes on the inner intervals; the outer
two never terminated) survive as corroboration only.

**The census is not a function of the three determinants alone.**  The
sign-chart bridge transfers verbatim (`trap ⟺ Wtau·Wpot < 0`,
`@positive ⟺ Wtau·Wwgt > 0`), but the noncentered census carries
`separate:` rows: each locate-then-certify list has a separated half whose
returned families are individually Krawczyk-certified with interval Schur
types. The exhaustive branch-and-bound closes that list at the twenty-two
named teachers. Away from those teachers, an unseeded family is not excluded
unless a regional atlas covers the point; `CERTIFICATE.md` sections 4 and 9
separate those scopes.

## The antipodal band, and the locator fix

The 1868 fragments the pre-fix map drew near `beta = pi` were a failure of
`separatedScan`: a uniform seed lattice stepping over the near-coincident
families (gaps 0.02–0.07) that carry the band's census, an `|h(D)| > 5e-3`
guard blind to the antipodal stratum `D = pi`, and a root-merging miss
within `3e-3` of `beta = pi`.  The certificate's replacement — weights from
the radial rows, whose determinant `pi^2 − phiJ(D)^2` vanishes only at
student coincidence, with a geometric gap grid — is certified per family
here and, since commit `2201dd0`, is the shipped locator
(`separatedScanJ`), formula for formula. The post-fix classifier agrees with
every accepted returned list in that band scan. The scan observes three
certified returned row-list patterns depending on `y` alone across its sampled
`beta` values; it is neither an exhaustive census there nor a proof of
constancy on the whole band.

## The sixth observed selected-row list

Hugging the strata `y = 0` and `|y| = 1` in the mixed sector, the locator
loses one of the proposed label's two `separate:trap@mixed` families. At
`(beta, y) = (0.15, 0.004)` and its exact mirror
`(0.15, 0.996)`, every located family in the resulting row list is
Krawczyk-certified and the coincident chart is complete. There is no exclusion
cover for the separated complement, so this is not an exact census.
`wall_pos.json` gives a two-sided locate-then-certify transition pair only
at `beta = 0.15`; its other eight sampled gaps, extending through
`beta = 3.00`, give one-sided first-success heights. No global collar
stratum, nine-face count, or wall continuation is certified
(`CERTIFICATE.md` 6.3/9.3).

## The certificate format, in full

A completeness certificate for one teacher is a **bisection tree** over the
`(th1, D)` domain, emitted by `face_bb_signlaw.py` and replayed by
`ivcert_check.py`.  The search that produced it — heap, priorities, ordering —
is prospecting and appears nowhere in the certificate.  This section specifies
the encoding completely: a referee should be able to write the reader from it.

The **spec** (`beta`, `y`, optional `beta_hi`/`y_hi`, `delta`, `seam`,
`dmax_pad`, `atom_eps`, `cell`, `minw`, `n_root_cells`) determines the root
cells, and the checker REGENERATES them from it; no box coordinate appears in
the file.  The **split rule** is determined by the box — bisect the WIDER side
(the `s`-side when `s1 − s0 ≥ d1 − d0`), child `0` is the lower half — so a
reader that knows the root cells reconstructs every box exactly, and no float
comparison is ever needed to identify one.

**`delta` is a validity domain, not a tuning knob.**  The certificate is a
completeness statement about gaps `D ≥ delta`, so at a named teacher `delta`
must sit below its smallest family gap, or the certificate is complete about a domain
that omits real families.  Read the `delta` of a certificate before comparing
its count with a census.

What the certified data says about the actual floors: `delta = 0.02` is below
every gap at the F2, F4, F5 and F6 witnesses, and it is ALSO below every gap
at the F3 and F7 witnesses — both certify their three separated families at
`D = 0.0315, 0.0597, 0.0961` (identical at the two, as the `y → −1−y`
teacher-swap mirror requires).  So the `delta = 0.005` those two runs used was
not needed at those witnesses, and a rerun at `0.02` would cover a smaller
domain and cost less.

**RE-CORRECTED, and this is the substantive version.**  An earlier note here
claimed that `delta = 0.02` excludes real families at `beta ~ 3.1`; last round
I retracted that on the strength of the F3 and F7 certificates, whose smallest
gap is `0.0315`.  Checking ONE `beta` and generalising was the mistake — the
retraction was wrong and the original claim was right in substance.

Tracking the three certified families of the F3 column by Newton continuation
as `beta -> pi`:

| `beta` | 3.08 | 3.10 | 3.102 | 3.104 | 3.12 | 3.13 | 3.14 |
|---|---|---|---|---|---|---|---|
| smallest gap | 0.0315 | 0.0212 | **0.0202** | **0.0192** | 0.0110 | 0.0059 | 0.00081 |

**The gaps vanish LINEARLY at the degenerate column**: measured
`D_min / (pi - beta) = 0.5082` and constant to four digits over three decades,
down to `pi - beta = 1e-6`.  All three families of the column collapse into the
coincident stratum as `beta -> pi`, which is a real feature of the map and not
a numerical artefact — `beta = pi` is exactly where the certificate's own
analysis says the separated stratum degenerates.

So the gap floor is not a tuning knob at all but a **named boundary of the
completeness statement**, with an explicit validity bound:

> a certificate with gap floor `delta` is complete only for
> `beta < pi - 2*delta` (measured crossings: `delta = 0.02` at
> `beta = 3.1024`, `delta = 0.005` at `beta = 3.1318`; the law predicts
> `3.1016` and `3.1316`).

**This boundary is now stated in the Lean face theorem, not only here.**
`censusAngleMapJ_familyCount_eq_witness_on_censusFace`
(`CountConstancyJ/SquareCountJ.lean`) takes the face as a rectangle with
explicit bounds and DERIVES from them the two side conditions the count law
needs, so every bound in its statement is used and none is decoration.  The
`beta`-upper bound is where the excluded neighbourhood of the degenerate
column lives.  The truncated face is still preconnected — excluding the column
only lowers the upper bound and the region stays a product of closed
intervals, so `isPreconnected_censusBoxJ` applies verbatim; the proof goes
through it, so that is checked in the proof term rather than asserted.

**No fixed `delta` covers a face up to `beta = pi`**, so the F3 and F7 face
statements will need the degenerate column excluded by a named neighbourhood
rather than by a smaller floor.  No existing certificate is invalidated: every
witness sits well inside its own bound (`beta = 3.08` against `3.1316` for
`delta = 0.005`; the `delta = 0.02` batch runs to `beta = 3.0` against
`3.1016`).

## Why the fold-freeness sweep does not scale to a face, and it is dimensional

The eight-worker authorisation was not spent, and this is the measurement that
says it should not be.  Descending one cell to a fixed depth under each MAP
form and counting boxes still alive there — the quantity that governs a
bisection sweep's total size, since a box alive at depth `k` spawns the whole
subtree below it:

| depth | plain row form | plain + mean-value | ratio |
|---|---|---|---|
| 8 | 178 alive | 108 alive | 1.6× |
| 12 | 1 448 | 412 | 3.5× |
| 16 | 11 824 | 3 594 | 3.3× |

Two readings, and the second is the one that matters.

**The mean-value MAP test is worth having but is a constant.**  The ratio
saturates near `3.3×` rather than compounding, against a per-node cost about
`2.4×` higher, so the net is modest.

**Both arms grow at about `1.7` survivors per LEVEL, and the irreducible floor
is `1.41`.**  The zero set of the angle map is codimension two in the
four-dimensional region `(beta, y, s, D)`, so covering it costs four boxes per
halving of every side — `4^{1/4} = 1.41` per level — no matter how good the
enclosure is.  Both forms sit just above that floor.  So the sweep's cost is
set by GEOMETRY, not by enclosure quality, and no improvement to the MAP test
can change the exponent: box count scales as `h^{-2}` in the target side `h`,
which for the `~1e-3` the MAP test needs, from a `D`-range of `3`, is of order
`1e7` boxes for ONE cell of a teacher rectangle only `0.002 × 0.002` wide.

**The direction that could change the exponent is closed, and by a vanishing
quantity rather than a dependency.**  The determinant test is `y`-free, so it
lives in three dimensions where the same codimension-two set costs `h^{-1}`
instead of `h^{-2}` — covering a curve rather than a surface.  Collecting that
requires the determinant test to fire in the small-gap regime, and it cannot.
Measured across the certified families, sorted by gap:

| `D` | `|det T|` | its slope | `|det|/|det′|` | levels needed |
|---|---|---|---|---|
| 0.359 | 6.4e-06 | 1.2e+00 | 5.4e-06 | 17 |
| 0.202 | 1.2e-13 | 1.6e-02 | 7.9e-12 | 37 |
| 0.096 | 1.6e-26 | 3.2e-08 | 5.2e-19 | 61 |
| 0.032 | **1.6e-36** | 1.7e-11 | **9.6e-26** | **83** |

As the gap shrinks by a factor 11, `|det T|` falls **31 orders** while its
derivative falls only **11**, so the ratio `|det|/|det′|` — the distance to the
nearest zero of the determinant — falls twenty orders, to `1e-25`.  A sign test
needs boxes smaller than that distance, which is 83 bisection levels.

That is a VANISHING quantity, not a dependency slope: no enclosure form, of any
order, can sign a quantity whose own zero sits `1e-25` away while its
derivative is `1e-11`.  The mean-value form had already removed the first-order
dependency here — it bought 20× to 1000× on exactly these families — and the
requirement stayed below `1e-7`.  What is left is the separated stratum's Schur
structure degenerating as the two students merge, the same degeneracy as at the
exact fit in graded form.

> **RETRACTED 2026-08-10, and the correction names a different repair.**  The
> two paragraphs above are wrong about the MECHANISM, though not about the
> verdict.  `|det|/|det′|` as tabulated is not the distance to a zero of the
> determinant.  Measured independently — central differences of the literal
> determinant at 160 bits, deliberately not through the mean-value form's
> forward-mode AD path, so the two computations share no code
> (`row7a_band.py`) — the linearised distance
> `|det| / (|∂_s| + |∂_D| + |∂_β|)` at those same families is `1.0e−3`,
> `9.9e−4`, `9.4e−4` at `D = 0.0315, 0.0597, 0.0961`, not `1e−25`.  It does not
> collapse with the gap, because `|det|` falls 31 orders across that range and
> its derivative falls with it: the determinant is a degree-four form in a
> kernel vector that itself goes to zero as the students merge, so numerator and
> derivative scale together.  **The determinant is not near a zero at the
> small-gap families.  The enclosure cannot see the sign, which is a different
> claim with a different repair** — stage by stage (`row7a_stagewidth.py`) every
> intermediate is exact to `1e−36` at a THIN point, `S0` and `N0` improve
> LINEARLY with the box side as a dependency should, and only the assembled
> determinant stalls, at a relative width near `40` below a side of `1e−8`.
>
> How much of the enclosure is dependency, measured against the determinant's
> sampled TRUE range over the same box (`row7a_stagewidth.py`, headroom
> section):
>
> | family | `D` | box side | `|det|` | true range | mean-value width | factor |
> |---|---|---|---|---|---|---|
> | F3 `(3.08,−0.08)` | 0.0315 | `1e−4` | `1.24e−34` | `1.22e−35` | `3.73e−14` | **`3.0e+21`** |
> | F3 `(3.08,−0.08)` | 0.0597 | `1e−4` | `2.39e−28` | `2.42e−29` | `4.63e−12` | `1.9e+17` |
> | F6 `(0.75,−0.98)` | 1.4530 | `1e−4` | `1.42e+04` | `3.29e+02` | `4.58e+05` | `1.4e+03` |
>
> Over all 55 non-exact-fit certified families at a fixed box side of `1e−4`
> (`row7a_depfactor.log`), that factor is a clean function of the STUDENT GAP
> and of nothing else visible:
>
> | `D` band | `0.00–0.05` | `0.05–0.15` | `0.15–0.30` | `0.30–0.50` | `0.50–1.00` | `1.00–2.00` | `2.00–4.00` |
> |---|---|---|---|---|---|---|---|
> | families | 2 | 4 | 3 | 8 | 5 | 15 | 18 |
> | dependency factor | `3e21`–`7e21` | `2e14`–`3e18` | `7e05`–`8e08` | `1e02`–`2e05` | `9e00`–`4e02` | `2e00`–`1e03` | `1e00`–`3e02` |
>
> — a log-log slope of `−9.17` across four decades of gap (residual scatter
> `2.3` decades).  **So the `D^9.9` law this file recorded for the reachability
> threshold is not contradicted: the required box side really does collapse like
> a high power of the gap.  What is corrected is where that power lives** — in
> the enclosure's dependency factor, not in the determinant's distance to a
> zero, which is flat at about `1e−3` across the same range.
>
> At the smallest-gap family the determinant varies by **10 %** of its value
> across a box of side `1e−4` — exactly what a distance-to-a-zero of `1e−3`
> predicts, which is an independent confirmation of that number by a method
> sharing no code with the central differences — while the enclosure is `3e21`
> times wider than that variation.  An exact-range evaluation would sign the box
> at a side of `1e−4` with a factor of ten to spare, at EVERY certified family.
> The dependency factor is what blows up as the students merge (`1.4e+03` at
> `D = 1.45`, `3.0e+21` at `D = 0.03`), not the determinant's distance to a
> zero.  So
> the repair class is a stage-wise centered (Taylor-model) arithmetic, in which
> each intermediate carries a thin centre and an interval gradient; the
> output-level mean-value form already in place is exactly the thing whose own
> gradient enclosure stalls.  What is NOT retracted: the determinant test still
> never fires in the small-gap band — measured inside the sweep itself, zero
> SCHUR kills over the whole gap band `[0.02, 3.19]` against 9 896 over
> `[1.00, 3.19]` on the same rectangle with the same code — and the covering
> cost is dimensional either way, so this correction does not reopen the route.
> It changes what a successor would have to build to move the boundary, and it
> removes an impossibility claim that was never established.

So the three-dimensional route cannot be collected where it is needed, and the
fold-freeness sweep has no remaining path to a face by this design.  What that
leaves is an obligation the sweep was only ever one way of discharging: the
Lean face theorem still takes fold-freeness as an explicit hypothesis, and a
different discharge — per-family rather than by covering, or the mass-carrying
route that already covers the exact fit — is a different problem from the one
attacked here.

**The superseded reading, kept because it was published:** an earlier version
of this passage named the `y`-free determinant test as the direction that could
change the exponent.  It is the right direction and the wall is not in it; the
wall is that the quantity it tests vanishes exactly where the sweep needs it.

## Pricing the three routes to fold-freeness

The Lean face theorem takes fold-freeness as an explicit hypothesis; three ways
of discharging it have now been priced, and they differ by five orders of
magnitude.

| route | what it establishes | measured price per face |
|---|---|---|
| covering sweep over `(beta,y,s,D)` | `hfold` at EVERY zero, directly | `~1e7` boxes for one cell of a `0.002 x 0.002` rectangle — years |
| parametrized enumeration over teacher boxes | the complete count on each box | boxes capped at width `0.008` by solution drift ⇒ `~39 000` full enumerations — years |
| **mass-carrying continuation** | `hfold` at the TRACKED zeros, with the Schur type certified at each | **~2–4 minutes** — now BUILT and RUN, see the next section for the artifact and the wall clock |

The third is measured, not estimated.  Carrying the certified F4 family by
Newton prediction plus Krawczyk certification, and reading the interval Schur
type at every step:

* in `beta`: steps of `0.05` certify, at `46–50 ms` per step, type `saddle`
  throughout;
* in `y`: steps of `0.05` certify along the stated track
  (`y = −0.40` to `−0.85`), at `47–54 ms` per step, type `saddle` at every
  certified step.

A step of `0.05` is **six times** the `0.008` ceiling that caps the
parametrized route, and for the reason that matters: a predictor re-certifies
at each new parameter, so no single box has to contain the solution for a whole
range of teachers.  Gridding a face at that step — `beta` over `(0,pi)`, `y`
over `~0.76`, about `1 000` points per family — is `~2 000` certified steps for
F4 and `~4 000` for a three-family face.

**What continuation does NOT give, stated plainly.**  It certifies
fold-freeness along the sheets it tracks.  The schema's `hfold` quantifies over
ALL zeros, and continuation cannot by itself exclude a zero that is not a
continuation of a witness family, nor the birth of a new pair.  There is a
bootstrap — if `hfold` holds at every zero the count is constant, so every zero
IS a continuation — but it is circular until a degenerate zero is excluded
independently.

**That residual is a smaller problem than the sweep, and by an exponent.**  A
degenerate zero satisfies three equations in four unknowns (`map = 0` twice,
`det = 0` once), so it is a CURVE in the region, and covering a curve costs
`h^{-1}` where the sweep's surface cost `h^{-2}`.

**But the two boundaries are NOT the same one, and that is measured.**  The
determinant test needs boxes smaller than `|det|/|det′|`; across all 50
certified families that scale follows `|det|/|det′| ~ D^{9.9}`, and
reachability at a box side of `1e-8` ends at

> `D* ≈ 0.40`,  against certified-domain floors of `delta = 0.02` and `0.005`
> — a factor of **20** and **80**.

So there is a wide band, `delta ≤ D ≲ 0.4`, INSIDE the certified domain where
no determinant-based test can fire, and the band is not vacuous: **eight of the
certified families live in it**, at `D = 0.0315, 0.0597, 0.0961` (F3 and F7,
twice each) and `D ≈ 0.20` (F2 and F6).  Every one of the four faces those
belong to has a family the determinant route cannot serve.

**The honest statement is therefore the conditional one.**  Fold-freeness is
dischargeable by the determinant where the families are well separated
(`D ≳ 0.4`), and below that it remains a genuine hypothesis of the Lean face
theorem rather than something this machinery can supply.  That is not a defect
of the certificates — they are complete about `D ≥ delta` and say so — it is
the boundary of what a determinant sign test can see, and it sits an order of
magnitude above the domain's own floor rather than coinciding with it.

## The mass-carrying continuation, built and measured

The third route of the table above is now an ARTIFACT, and every figure in
this section is a wall-clock measurement of it rather than a projection from a
step cost.  `cont_cert.py` is the engine and the emitter, `cont_check.py` the
replay checker, `cont_gate.py` the falsifier for the one new quantity.  Read
the scope paragraph at the end before quoting any of it: **this does not
discharge the Lean face theorem's `hfold`, and its own checker says so on
every successful run.**

### What one step is

A step carries one tracked zero of the census angle map from a teacher at
which it is already certified to a neighbouring teacher.

1. **PREDICT** — a damped Newton corrector for the row form `(E0, E1)` in the
   `(s, D)` chart, at the working precision on thin intervals, from the
   midpoint of the enclosure certified at the previous teacher.  A guess and
   nothing else: it enters no verdict, and a bad guess can only make a
   certification fail.
2. **CERTIFY** — the same `2 × 2` Krawczyk test the enumeration uses
   (`census_cert._krawczyk_sD`) on the box of radius `r` about the predicted
   point, `r` halved from `r0` until `K(Z) ⊆ int(Z)` — exactly one zero of the
   map in `Z`, enclosed — or `r < rmin`, which is an UNDECIDED verdict.
3. **READ THE TYPE** — on the Krawczyk-refined enclosure: the interval Schur
   type (`census_cert.schur_label`, the `widgets.js` label in interval
   arithmetic), and the sign of the ANGLE-MAP JACOBIAN DETERMINANT.

### The new quantity, and its gate

The determinant the Lean schema consumes is `separatedAngleJacDetJ`
(`SeparatedIndex/AngleJacJ.lean`), the Jacobian determinant of the
angle-equation map in `(θ₀, θ₁)`.  In the `(s, D)` chart the row form's own
derivative (`census_cert.sep_res_rowform_jac`) gives it with no new
arithmetic,

```
   det_theta(E0,E1)  =  dE0/dD · dE1/ds  −  dE0/ds · dE1/dD
```

because `d/dθ₀ = d/dD` and `d/dθ₁ = d/ds − d/dD`.  The tree's quantity is
`separatedAngleJacDetJ = −det_theta(E0,E1)`, and the sign is derived rather
than chosen: the tree's second equation carries `couplingHJ (θ₁ − θ₀)` and
`couplingHJ` is odd, so it is the NEGATIVE of this directory's `E1`, and one
negated row flips the determinant.  A SIGN-DEFINITE enclosure of `det_theta`
over a Krawczyk-certified solution box is therefore `separatedAngleJacDetJ ≠ 0`
at the enclosed zero — fold-freeness at that zero, with no detour through
either the mass-free type boundary or the Schur block.

`cont_gate.py` is the falsifier: three checks over 200 random teachers and
student pairs at 300 bits, all evaluated at the SAME exact `(θ₀, θ₁)` (rounding
`θ₀ = s + D` through a floating-point sum alone shows up at `1e-16` and would
mask a real defect).

| check | worst over 200 samples |
|---|---|
| `det_theta` against central differences of the row form in `(θ₀, θ₁)` | `3.3e-49` |
| `det_theta + separatedAngleJacDetJ`, from a transcription of the tree's four closed forms | `3.4e-88` |
| those closed forms against central differences of their own equations | `1.0e-49` |

An invalid derivative enclosure is the failure this certificate has already
paid for once (the 2026-08-08 Krawczyk sign fix), so the quantity is gated
before any leaf claims fold-freeness with it.

### Why a bisection tree, and over what

A continuation looks like a flat list of steps, and packing a list into a
bisection tree would be packaging rather than mathematics.  It is not a list.
When a step is refused the walk HALVES the parameter interval and goes through
the midpoint, so the steps of one line ARE the leaves of a bisection tree over
that line's parameter interval, and **preorder traversal is exactly walk
order**.  A node is a parameter interval `[u₀, u₁]` in the line's normalised
coordinate: `1` = internal, the step was refused and the walk goes
`u₀ → m → u₁` through the midpoint `m = (u₀+u₁)/2`; `0` = leaf, the walk
stepped straight to `u₁` and the three following bits are the verdict there.
The midpoint rule is determined by the node, so no coordinate appears in the
file and the checker regenerates every parameter from the spec — the same
discipline as `face_bb_signlaw.py`.  Root cells are the GRID intervals, so
every grid teacher is the right endpoint of some leaf and carries a verdict;
bisection adds certified teachers between grid points and never removes one.
Bits are packed MSB-first and base64-encoded as `bits_b64`, with `n_nodes`
giving the meaningful bit count, exactly as in `ivcert_pack.py`.

The lines, in bitstream order, for each tracked family in turn: the ANCHOR (a
degenerate leaf at the first grid teacher, from that family's declared float
seed), the SPINE (a walk in `y` at the first `β`), then one TOOTH per `y`-grid
value (a walk in `β`, anchored at the spine's certified state there).  The
seeds are spec INPUT, like `delta` and `seam`: they are located by the float
locator, which is prospecting, and a wrong seed can only make a leaf fail.

| code | verdict | what the leaf claims |
|---|---|---|
| 0 | `CERT_SADDLE` | one zero, enclosed; `separatedAngleJacDetJ ≠ 0`; Schur type `saddle` |
| 1 | `CERT_SPURIOUS` | the same, Schur type `spurious` |
| 2 | `CERT_FOLDFREE` | the same, Schur type undecided |
| 3 | `FOLD_CANDIDATE` | one zero, but the Jacobian-determinant enclosure straddles zero |
| 4 | `UNDECIDED_JUMP` | certified, but further than `maxjump` from the previous enclosure |
| 6 | `UNDECIDED_MINRAD` | Krawczyk never certified down to `rmin` |
| 7 | `UNDECIDED_OTHER` | any other reason |

Only 0, 1 and 2 are accepting, so an undecided leaf fails the replay exactly as
in `ivcert_check.py`.  Codes 3 and 4 are kept distinct from 6 and 7 because
they say something: 3 is a fold candidate, 4 is a walk that left its sheet.

### The first F4 run was WRONG, and its own disjointness check caught it

Recorded because it is the failure mode of the method, not an accident of this
run.  The first F4 certificate was emitted with an undamped corrector and no
sheet guard.  Every one of its 1 506 leaves was valid and re-derived
identically on replay — and the exact-fit walk had stepped onto the separated
family's sheet at `(β, y) = (0.20, −0.65)`: its "family 1" box was family 0's
box with the labels swapped and shifted by `2π`.  The certificate was
internally sound and its READING as two tracked families was false.

**In a continuation the per-step verdict cannot detect a sheet jump** — every
box it certifies is a genuine box around a genuine zero — **so the tracking
has to be checked by a separate global invariant.**  Two now exist, and both
are part of the deterministic step, hence enforced by the checker:
end-to-end pairwise disjointness of the families AS UNORDERED PAIRS (the swap
`(s,D) → (s+D, 2π−D)` applied before the comparison), and a per-step sheet
guard.

### The three rules that followed, each with the measurement that forced it

* **A trust radius plus backtracking on `|E0| + |E1|`.**  The undamped step is
  `J⁻¹E`, and `det J` at the exact fit collapses as the fit merges into the
  coincidence stratum at the degenerate column `β = 0`.  Measured at
  `y = −0.65`, `separatedAngleJacDetJ` at `θ = (β, 0)`:

  | `β` | 0.05 | 0.10 | 0.15 | 0.20 | 0.30 | 0.60 | 1.00 | 2.20 |
  |---|---|---|---|---|---|---|---|---|
  | `det` | `8.8e-7` | `5.5e-5` | `6.0e-4` | `3.3e-3` | `3.4e-2` | `1.63` | `20.8` | `293` |
  | slope `d log|det| / d log β` | | 5.96 | 5.92 | 5.88 | 5.80 | 5.49 | 4.98 | 2.66 |

  so `det ~ β^5.9` at the column, and a full Newton step from the exact fit at
  `β = 0.15` toward `β = 0.20` has length `3.2` and lands on another family.
  **This is a conditioning wall of the PREDICTOR, not of the sign test**, and
  the lane's own diagnostic separates them: the box side a sign test needs is
  `|det| / |det′| ≈ β/6`, LINEAR in `β`, against an enclosure width of `~1e-40`
  at a refined solution box — the determinant is signed with enormous margin at
  every `β` in the table.  A vanishing quantity is one kind of wall for a SIGN
  test and a different, milder one for a SOLVER; this lane had only ever
  measured the first.
* **A gap floor `delta` in the corrector.**  The row form `E = massDet(D)·G`
  vanishes IDENTICALLY on `D = 0` (`massDet(0) = 0` and `h(0) = 0`), so the
  coincidence line is a SPURIOUS zero manifold of the row form and an undamped
  corrector is attracted to it — measured: the exact-fit step
  `β = 0.10 → 0.15` collapses to `D ≈ 1e-10`.  `delta` here is the same
  validity-domain object as the enumeration's gap floor, not a tuning knob.
* **The sheet guard `maxjump`.**  A certified box further than `maxjump` from
  the previous one is a step that left its sheet; the verdict is
  `UNDECIDED_JUMP`, which makes the walk bisect the parameter interval and
  retry on half the step.  The box it refuses is a perfectly good enclosure of
  a perfectly good zero; what fails is the reading of the sequence as ONE
  tracked family, and that reading is the whole point.

### F4, end to end

`cont_F4.json`, 742 bytes packed, replayed by `cont_check.py`
(`emit_cont_F4.log`, `replay_cont_F4.log`).  The rectangle is
`β ∈ [0.10, 2.20]` and `y ∈ [−0.90, −0.10]`, both at the priced step `0.05`.
The proposed F4 label is `0 < β < β*₁ = 2.2257` outside the torque lens and
`−1−w_a(β) < y < w_a(β)`; a located family count of 1 (against 3 in F2 above
and F6 below) is observed for `y ∈ [−0.91, −0.09]` at every `β` probed, so this
grid samples the proposed label less a collar of about `0.01` at each proposed
`y` boundary and `0.1` / `0.026` at its two `β` ends. Those candidate
boundaries are sampled rather than enclosed, so the rectangular grid is not a
certified global cell.

| | measured |
|---|---|
| grid | `43 × 17 = 731` teachers |
| tracked families | 2 — the separated family of F4, and the exact fit `θ = (β, 0)` |
| certified steps | **1 479** (731 per family plus 17 bisection midpoints) |
| parameter bisections | 17, all at the small-`β` end |
| undecided leaves | **0** |
| verdicts | `CERT_SADDLE` 731, `CERT_SPURIOUS` 748 |
| disjointness | certified at every grid teacher, worst separation `1.541` |
| worst step displacement | `0.1622`, against the guard `0.35` |
| unordered gaps over the grid | family 0 in `[1.5999, π]`, family 1 in `[0.10, 2.20]`, both above `delta = 0.02` |
| **emit** | **122.8 s** (83.0 ms per certified step) |
| **replay** | **121.4 s** (82.1 ms per re-derived step) |

Two minutes, single core, for a rectangle spanning much of the proposed F4
label — against
`~1e7` boxes for one `0.002 × 0.002` cell by the covering sweep and
`~39 000` full enumerations by the parametrized route.  The replay costs what
the emit costs, because the checker re-derives every verdict rather than
reading it; the certificate's economy is in the 742 bytes, not in the replay.

The separated family is `saddle` at all 731 teachers and the exact fit is
`spurious` at all of them, which is the same interval Schur reading that row 7b
records at the 22 replayed witnesses — here at 731 grid teachers in one rectangle
instead.

### F6, four families

`cont_F6.json` (`emit_cont_F6.log`, `replay_cont_F6.log`).  F6 is THIN: a
located family count of 3 holds only for `y` between about `−0.915` and `−1`,
so the `y` step here is `0.02` and not `0.05`.  Rectangle `β ∈ [0.30, 2.20]`
step `0.05`, `y ∈ [−0.98, −0.94]` step `0.02`; four tracked families — F6's
three separated families and the exact fit.

| | measured |
|---|---|
| grid | `39 × 3 = 117` teachers |
| tracked families | 4 — F6's three separated families and the exact fit |
| certified steps | **578** (468 grid leaves plus 110 bisection midpoints) |
| parameter bisections | **110** — F6's families run close together (worst separation `0.2818`), so the sheet guard fires often and the walk halves its step |
| undecided leaves | **0** |
| verdicts | `CERT_SADDLE` 238, `CERT_SPURIOUS` 340 |
| disjointness | certified at every grid teacher, worst separation `0.2818` |
| worst step displacement | `0.2359`, against the guard `0.35` |
| unordered gaps over the grid | `[0.113, 0.659]`, `[0.723, 1.605]`, `[1.500, 3.064]`, `[0.300, 2.200]`, all above `delta = 0.02` |
| **emit** | **73.3 s** (126.8 ms per certified step) |
| **replay** | **56.6 s** (98.0 ms per re-derived step) |

The per-step cost is higher than F4's because a refused step costs a
certification that no leaf records: with 110 of 578 leaves reached after a
bisection, the emitter pays for the attempts and the replay does not, which is
why the replay is the cheaper of the two here and was the equal of the emit on
F4.

### F1, the seventh proposed label — five tracked families on grids spanning `0.30 ≤ beta ≤ 3.00`

F1 was the one proposed label without a continuation certificate, and the reason was structural: the
covering sweep is what had failed there, and the enumeration's order-only
prefilter is weakest in the mixed sector (it discards 39 % of the area at the
F1 witness `(2.60, 0.20)` against 81 % in the positive one).  **The
continuation does not use the sweep at all**, so that obstruction does not
transfer, and F1 now has the LARGEST continuation certificate
here: `0.30 ≤ beta ≤ 3.00` and `0.10 ≤ y ≤ 0.90` at step `0.05` in both, five
tracked families — F1's four separated families and the exact fit, more than
any other proposed label — over `55 × 17 = 935` teachers.

It is emitted as two rectangles that share the column `beta = 2.60`, for the
reason in the design-rule section below and for no other: F1's family
separation falls by a factor of seven across the rectangle, and one sheet guard
cannot serve both ends.

| | `cont_F1a.json` | `cont_F1b.json` |
|---|---|---|
| rectangle | `beta ∈ [0.30, 2.60]`, `y ∈ [0.10, 0.90]` | `beta ∈ [2.60, 3.00]`, `y ∈ [0.10, 0.90]` |
| grid | `47 × 17 = 799` teachers | `9 × 17 = 153` teachers |
| tracked families | 5 | 5 |
| certified steps | **9 643** | **10 219** |
| parameter bisections | 5 648 | 9 454 |
| undecided leaves | **0** | **0** |
| verdicts | `CERT_SADDLE` 4 880, `CERT_SPURIOUS` 4 763 | `CERT_SADDLE` 5 399, `CERT_SPURIOUS` 4 820 |
| disjointness | certified at every grid teacher, worst separation `0.1152` | certified at every grid teacher, worst separation `0.02386` |
| sheet guard / worst displacement | `0.04` / `0.03999` | `0.008` / `0.008` |
| gap floor `delta` | `0.005` | `0.005` |
| packed size | 5 528 bytes | 6 292 bytes |
| **emit** | **1 292.3 s** | **2 142.8 s** |
| **replay** | **895.7 s**, VALID AND COMPLETE | **1 077.5 s**, VALID AND COMPLETE |

**Where the rectangle's four walls come from — all four are a vanishing, none
is a preference.**  This is the `0d` test applied to the certificate's own
thresholds: not "is the floor small enough at the witness I tested", but "does
the thresholded quantity have a ZERO inside the region I am claiming".  It
does, at all four walls, and that is why they are named boundaries of the
claim rather than parameters in a script.

* `y → 0` and `y → 1`, the two teacher-mass strata.  F1's smallest family gap
  vanishes LINEARLY at each: `D_min / y` is `0.115` at `beta = 0.30`, `0.26`
  at `0.50`, `0.53` at `0.75`, `0.87` at `1.00`, measured down to `y = 0.008`
  (`work/contF1/scan1.log`; the `locator_failed` rows of the older
  `wall_pos.json` are the same family failing to certify as its gap closes).
  **No positive gap floor covers a rectangle that reaches either stratum.**
  The collar is `0.10` wide because at `y = 0.10` the worst gap over the whole
  `beta` range is still `0.0124`, i.e. `2.5×` the floor used.
* `beta → pi`, the degenerate column.  The same statement in the other
  coordinate: `D_min / (pi − beta) = 0.66` at `y = 0.5`, flat over
  `beta = 3.05` and `3.10`, and the four families crowd together with it
  (pairwise separation `0.0182` at `beta = 3.00` against `0.132` at `0.30`).
  At `beta = 3.10` the four gaps are `0.0056 / 0.0095 / 0.0397 / 0.0573` at
  `y = 0.05` and the closest pair is `0.0028` apart.
* `beta → 0`.  As on F4, the exact fit merges into the coincidence stratum and
  the wall is the conditioning of the walk, not the reach of the enclosure.

The Schur reading is the same at all 935 grid teachers. It matches the exact
census at the exhaustively enumerated teacher `(0.75,0.40)` and the
Krawczyk-certified returned lists at the other three F1 witnesses: two
`saddle` separated families, two `spurious` ones, and a `spurious` exact fit.

### F1's per-teacher ENUMERATION certificate, which is a different claim

The continuation above tracks five families at the grid teachers spanning the
proposed-label rectangle. Separately, and for the first time at a teacher assigned this
proposed label, the sign-law-pruned enumeration ran to
completion at F1's representative witness `(0.75, 0.40)`:

| | measured |
|---|---|
| domain area | `19.9277` |
| steps | **509 533** |
| leaves | **254 893** (`ivcert_C_F1_0p75_0p4.bits.json`, 207.9 KB packed from 9 698 KB, `47×`) |
| undecided boxes | **0**, undecided area **0** |
| classes | 6 raw boxes → 5 classes, disjointness certified |
| separated families | **4**, plus 1 exact-fit class |
| wall clock, emit | **2 370.1 s** |

The four families are `D = 1.826063` (`saddle`), `0.580570` (`spurious`),
`2.424925` (`saddle`), `1.072995` (`spurious`), which is the census string
recorded for F1 at that teacher — now with the domain closed rather than
sampled.

**These two certificates say different things and neither implies the other.**
The enumeration says: at THIS teacher, the domain contains exactly these
classes and nothing else — completeness at a point.  The continuation says: at
935 teachers, these five sheets each carry exactly one zero with a nonzero
Jacobian determinant — fold-freeness along tracked sheets at those grid
teachers.
`(0.75, 0.40)` is a grid teacher of `cont_F1a.json`, so at that ONE teacher
the two coincide and `hfold` genuinely holds at every zero of the domain.  It
does not propagate, for the reason the closing section gives.

### F5 has differing exact point counts and measured fold evidence

The count change among teachers assigned proposed F5 is bracketed at
`beta = 3.0` by a pair of
complete enumeration certificates (three families at `y = −0.42`, one at
`y = −0.40`). The continuation records Jacobian determinants along selected
families and supplies evidence for a fold mechanism, but it does not certify a
degenerate intermediate zero.

Seed the three families that the locator finds on the equal-mass line
`y = −1/2` and walk them upward in `y` by certified Krawczyk steps of `0.002`
(`work/contF1/f5_foldscan.py`, whose run is `f5_foldscan.log`).  At every
`beta` tried, two of the three have the same last accepted `y` while the third
walks on to `y = −0.30` untroubled:

| `beta` | last certified `y` | `D` of the two tracked families | `det_theta` there |
|---|---|---|---|
| 2.30 | `−0.4980` | `3.0722` / `2.8027` | `−19.09` / `+19.81` |
| 2.60 | `−0.4580` | `2.5376` / `2.4975` | `−3.072` / `+2.973` |
| 2.90 | `−0.4240` | `2.3819` / `2.1120` | `−15.45` / `+10.71` |
| 3.00 | `−0.4160` | `2.2372` / `2.1506` | `−3.897` / `+3.424` |
| 3.10 | `−0.4140` | `2.3113` / `2.0100` | `−14.47` / `+8.956` |

The two tracked families have opposite determinant signs, and the recorded
determinant magnitudes fall along parts of the walk. This is consistent with a
saddle-node but does not certify one. At `beta = 2.60` the antipodal
branch's determinant runs `−117.5, −105.1, −86.2, −61.8, −23.2, −3.07` at
`y = −0.498, −0.490, −0.480, −0.470, −0.460, −0.458`, i.e. `det²` linear in
`y` under a numerical fit whose extrapolated zero is `y = −0.45796`, stable to
`1e−5` whether fitted on the last two steps or the last four. That fitted zero
is not enclosed. At `beta = 3.00`, the last certified height `y = −0.4160`
and the fitted zero `−0.41591` sit inside the
enumeration certificates' independent bracket `(−0.42, −0.40)` — two
compatible instruments with different evidence classes.

**The stopping height is instrument-dependent; the individual steps are not**,
so the instrument is stated with the table (`f5_foldscan.log`, regenerated
2026-08-10).  Steps of `0.002` in `y`, certify settings from `SPECS["F5b"]`
(160 bits, `r0 = 1e−3`, `rmin = 1e−11`, `trust = 0.04`), and the SHEET GUARD
`maxjump` DISABLED — deliberately, because `maxjump` is a tracking device
whose own design rule requires it to sit below the tracked families' measured
pairwise separation. With the guard off, every walk in the table stops on
`UNDECIDED_MINRAD`: the contraction test reaches its radius floor, which is an
algorithmic failure and not a degeneracy certificate. The two tracked walks
have the same last accepted height at all five betas. With `F5b`'s
`maxjump = 0.10` in place the walk stops one or two steps early at four of the
five betas on `UNDECIDED_JUMP` — a zero WAS certified there and only the sheet
reading was refused — the two branches part by a step at `beta = 3.00`
(`−0.4180` / `−0.4200`) and `3.10` (`−0.4160` / `−0.4180`), and at
`beta = 2.30` the FIRST step is refused, which is a statement about the guard
and not about the map.  Two numbers this file
carried before that run are withdrawn by it: the `det²` zero `−0.4585`, an
extrapolation from a walk that stopped at `−0.4600` and that sits BELOW a
height at which both families are certified; and the `beta = 3.00` wall
  `≈ −0.419`, which was the guarded instrument's stopping estimate. The per-step displacement
`|Δs| + |ΔD|` is printed at every step in place of the guard, so a sheet jump
stays visible: it never exceeds `0.209` on any walk here.

Three consequences, and the third is the one that matters for the Lean side.

* The pointwise count difference is real, not a locator artifact: complete
  enumerations give one and three families at nearby teachers. The determinant
  trends along the tracked sheets are evidence for a saddle fold, but the
  stopping heights do not enclose a global curve.
* The conditional `hfold` schema cannot simply be asserted on the
  undivided proposed F5 label. The differing point counts force at least one
  count-transport hypothesis to fail on a connecting region; the current
  evidence does not promote the proposed label to a certified face.
* The measured `|f| / |f′|` diagnostic returns its second answer. A large
  `|f′|` says the enclosure is at fault and a better form helps; `|f|`
  collapsing faster says the function's own zero is approaching and nothing
  helps. Here the sampled determinant trend suggests an approaching zero, but
  that interpretation remains numerical until an interval certificate encloses
  the degenerate zero. The two continuation rectangles therefore stay separate
  without being promoted to certified sides of a fold.

`cont_F5a.json` and `cont_F5b.json` track different parts of the
proposed F5 label. The split reflects the observed count change rather than
certifying a global subdivision.

### The design rule the proposed-label continuations produced

**The sheet guard must sit BELOW the family separation on the rectangle it
runs on, and the family separation is a measured property of the tracked
rectangle, not a
constant over a chosen rectangle.**

The guard `maxjump` refuses a certified step whose enclosure is further than
`maxjump` from the previous one, on the reading that the walk left its sheet;
a refused step halves the parameter interval and goes through the midpoint.
F4 and F6 ran at `0.35` and never needed more. The other five labels did, and
the price is recorded here rather than tuned quietly:

* **F3 passed at `maxjump = 0.05` and F7 failed at the same setting**, two
  tracked walks landing on ONE zero at `(2.97, −0.90)` with separation
  `−7.6e-44`. The two proposed-label grids are exact mirrors — `y → −1−y` composed with
  `theta → beta − theta` preserves the gap `D` identically — so F3's pass at
  that setting was luck, and reporting it as evidence would have been the
  error.  F7 failed again at `0.02`.  Both pass at `0.006`.
* The reason is a RATIO, not a threshold.  Over one grid step of
  `beta = 0.01` near the degenerate column the families move about as far as
  they are apart: `D` runs `0.0602, 0.0402, 0.0271` at `beta = 3.05, 3.08,
  3.10` while the closest pair there is `0.0126` apart.  No guard at or above
  the step displacement can separate them; below it, the bisection does the
  work and costs only wall clock (F3: 193 splits at `0.02`, 2 171 at `0.006`).
* **F1 forced the rule to be applied twice on one proposed label.** Its four families'
  worst pairwise separation runs `0.132` at `beta = 0.30`, `0.088` at `2.60`,
  `0.050` at `2.80` and `0.018` at `3.00` — a factor of seven across the
  rectangle — so a single guard is either useless at the crowded end or
  ruinous at the open one.  Measured, twice: at `0.15` the walk emitted 4 942
  leaves that all replay bit-exactly and still failed the disjointness check
  (families 2 and 3 on one zero at `(2.95, 0.30)`); at `0.05` it emitted
  9 497 and failed at `(2.95, 0.65)`.  F1's certificate is therefore two
  rectangles whose union is the original one, at guards `0.04` and `0.008`.

Two things worth separating, because they were confused once during this
pass.  A failed disjointness check does NOT mean a leaf is wrong: on every one
of those runs the replay re-derived all leaves bit-exactly and reported zero
failures and zero undecided.  What fails is the reading of a sequence of boxes
as ONE tracked family, and that reading is what the certificate is FOR.  And
the per-step verdict is structurally unable to see it — only the end-to-end
disjointness check can, which is the lesson F4 recorded on its first run and
which bit three more times here.

### Can the ORDER-ONLY exclusion decide the DETERMINANT? No, and the test is free

The sign law of `separatedNumJ` (`SeparatedCount/SignLawJ.lean`, axiom-clean)
licenses an ORDER-ONLY box exclusion — four endpoint comparisons of a student
offset against the student gap, no interval evaluation of anything — and that
prefilter is what made these enumerations tractable at all (it discards 15.89
of F4's 19.61 domain area, and 2.64 of F1's 19.93).  The obvious next thought
is to point the same prefilter at the ANGLE-MAP JACOBIAN DETERMINANT: if the
order data decided `sign(det)` it would exclude degenerate zeros over a region
without any reference to tracking, and the count bootstrap would stop being
circular.

**It does not; the F5 continuation supplies a counterexample.** The order data available to such a
test is the 4-vector

> `P = ( th0 vs D, th0 − beta vs D, −th1 vs D, beta − th1 vs D )`, each entry
> `<` or `>` mod `2π` — by the sign law exactly the sign of `separatedNumJ`
> at that argument.

**Read those four arguments off the Lean statements, not off memory.**  The
first version of this test used `th1` and `th1 − beta` for the last two; the
tree's `generalJAngleEq1_unitTeacher0/1` carry `−th1` and `beta − th1`, and
the wrong pair was caught only by re-reading the declarations.  The verdict is
unchanged and the corrected patterns are cleaner, which is luck, not a reason
to have skipped the check.

At the teacher `(beta, y) = (2.60, −0.480)` two separated critical
configurations on the F5 continuation are

| configuration | gap `D` | pattern `P` | `det_theta` |
|---|---|---|---|
| antipodal branch | `2.970075` | `><><` | `−86.2459` |
| mirror-pair branch | `2.056670` | `><><` | `+40.1637` |

— **the same order pattern and opposite determinant signs, at the same
teacher** (`work/contF1/detsign_test.py`; both are Krawczyk-certified zeros of
the census angle map, verdict `CERT_SADDLE`).  This is not one coincidence at
one teacher: the two branches carry the SAME pattern `><><` at every visited
teacher of the `0.005`-step walk, `y = −0.500` through `−0.475`, with
determinants of opposite sign at those visited teachers, while the third walk
continues to certify a family with the complementary pattern `<><>`. So the
order data is neither sufficient to fix the sign of `det` nor able to predict
which walk will stop certifying; a failed step makes no claim about the
underlying family after that teacher, and no
order-only exclusion built on it can deliver birth-freeness over a region.

Two things this does NOT say.  It says nothing against the sign law, which is
about `separatedNumJ` and is untouched; and it does not refute an interval
branch-and-bound on the determinant itself — that is the covering sweep, still
priced at `~1e7` boxes per `0.002 × 0.002` cell.  What it removes is the cheap
route to region completeness, and the reason it removes it is worth keeping:
the prefilter is cheap precisely because it reads only an ORDER, and a
determinant is not an order.

### Continuations on nine proposed-label rectangles

Every proposed label has a continuation rectangle, emitted by
`cont_cert.py` and replayed by `cont_check.py`, every one with **zero
undecided leaves** and certified pairwise disjointness of its tracked families
at the grid teachers. These are tracked-sheet certificates, not a cover of
global arrangement faces and not an exclusion of untracked families. Two
labels use two rectangles each — F1 because its
family separation varies by a factor of seven and one sheet guard cannot serve
both ends, and F5 because its exact point counts differ and the two
continuation rectangles certify different returned lists at their visited
teachers. This does not certify a
fold curve between the rectangles.

| proposed label | rectangle (`beta` × `y`, step) | teachers | tracked | certified steps | bisections | guard | emit | replay |
|---|---|---|---|---|---|---|---|---|
| F1a | `[0.30, 2.60]` × `[0.10, 0.90]`, `0.05` | `47 × 17 = 799` | 5 | 9 643 | 5 648 | `0.04` | 1 292.3 s | 895.7 s |
| F1b | `[2.60, 3.00]` × `[0.10, 0.90]`, `0.05` | `9 × 17 = 153` | 5 | 10 219 | 9 454 | `0.008` | 2 142.8 s | 1 077.5 s |
| F2 | `[0.30, 2.20]` `0.05` × `[−0.06, −0.01]` `0.01` | `39 × 6 = 234` | 4 | 939 | 3 | `0.35` | 87.7 s | 88.2 s |
| F3 | `[2.90, 3.10]` `0.01` × `[−0.10, −0.02]` `0.01` | `21 × 9 = 189` | 4 | 2 927 | 2 171 | `0.006` | 564.4 s | 324.0 s |
| F4 | `[0.10, 2.20]` × `[−0.90, −0.10]`, `0.05` | `43 × 17 = 731` | 2 | 1 479 | 17 | `0.35` | 122.8 s | 121.4 s |
| F5a | `[2.40, 3.00]` × `[−0.40, −0.15]`, `0.05` | `13 × 6 = 78` | 2 | 285 | 129 | `0.10` | 33.5 s | 22.5 s |
| F5b | `[2.60, 3.10]` `0.05` × `[−0.52, −0.48]` `0.005` | `11 × 9 = 99` | 4 | 914 | 518 | `0.10` | 108.8 s | 65.4 s |
| F6 | `[0.30, 2.20]` `0.05` × `[−0.98, −0.94]` `0.02` | `39 × 3 = 117` | 4 | 578 | 110 | `0.35` | 73.3 s | 56.6 s |
| F7 | `[2.90, 3.10]` `0.01` × `[−0.98, −0.90]` `0.01` | `21 × 9 = 189` | 4 | 3 294 | 2 538 | `0.006` | 663.3 s | 364.2 s |

**2 589 grid teachers, 30 278 certified steps, zero undecided leaves,
85 minutes of emit on one core** (2 572 DISTINCT teachers — F1's two
rectangles share the column `beta = 2.60`, and it is certified in both) — against `~1e7` boxes for one `0.002 × 0.002` cell by the
covering sweep and `~39 000` full enumerations by the parametrized route, both
of which are years.  Every wall-clock figure here is a measurement of the run
that produced the committed file, never a projection from a step cost; F1a,
F1b, F3, F5a, F5b, F7 and the F1 enumeration were emitted with between two and
five single-core Python processes running concurrently on a 24-core node, so
their per-step costs (`117–210 ms`) sit above F4's and F6's uncontended
`83–127 ms` and should be read as upper bounds.

A refused step costs a certification that no leaf records, so a continuation with many
bisections pays more to EMIT than to replay (F3: 564 s emit, 324 s replay);
one with few pays about the same for both (F4: 122.8 s and 121.4 s).  The
economy of these certificates is in the packed bitstring — 5 528 bytes for
F1a's 9 643 steps — not in the replay, which re-derives every verdict rather
than reading it.

The rectangles were chosen as strict sub-rectangles of the proposed labels.
The renderer draws candidate boundaries as curves through sampled windows,
but no such curve is globally enclosed, so no rectangle is itself a certified
global cell; each spec's comment names
which collar it drops and which measurement set guided the choice.

### Where the continuation and the enumeration MEET, which is at points

Not a coincidence and worth stating, because it is the only place the two
certificate kinds touch.  EIGHTEEN of the twenty-two enumerated teachers are
grid teachers of their proposed-label continuation rectangle — all of F1's, F2's,
F3's, F4's, F6's and F7's, and three of F5's.  The four that are not are
`(2.40, −0.50)` and the three point-count-bracket teachers `(3.00, −0.42)`,
`(3.00, −0.58)`, `(3.00, −0.60)`, which lie outside both F5 rectangles.  At those teachers the enumeration certificate says
the tracked list is COMPLETE — every point of the domain is excluded or inside
a Krawczyk box — so there the continuation's "the tracked zeros" and the
schema's "all zeros" coincide, and `hfold` genuinely holds at every zero.

**That does not propagate, and the temptation to say it does is exactly the
circular bootstrap.**  Completeness at a witness plus fold-freeness along the
tracked sheets gives the count at the neighbouring teacher ONLY IF no zero is
born in between, which is what fold-freeness at all zeros would have given us.
The two certificates meet at points; between the points they do not compose.

F5 is the standing reason this caution matters. Exact point counts differ
between teachers assigned the same proposed label, while the continuation
supplies only measured fold evidence. A bootstrap that treated the proposed
label as a certified face would therefore manufacture an unsupported face law.

### What this does NOT give, stated plainly

The `hfold` hypothesis of `censusAngleMapJ_zeroCount_eq_witness` and of
`censusAngleMapJ_familyCount_eq_witness_on_censusFace` reads

> for every teacher `p` of the region and every `x` of the compact student
> region with `censusAngleMapJ p x = 0`, `separatedAngleJacDetJ … ≠ 0`

— quantified over ALL zeros at ALL teachers.  A continuation certifies
fold-freeness **along the sheets it tracks, at the grid teachers only**.  It
cannot by itself exclude a zero that is not a continuation of a seeded family,
nor the birth of a new pair between two grid teachers.  There is a bootstrap —
if fold-freeness held at every zero the count would be constant, so every zero
would be a continuation of one at the witness — and it is CIRCULAR until
degenerate zeros are excluded independently.  Nothing here excludes them.

**The gap is named in the tree, not only here.**
`CountConstancyJ/TrackedFoldFreeJ.lean` writes `hfold` as the conjunction of
`CensusZerosContainedIn p K L` (every zero of the map in `K` at the teacher `p`
is one of the tracked ones) and `CensusJacDetNeZeroOn p L` (the Jacobian
determinant is nonzero at every tracked point), with
`censusAngleMapJ_foldFree_of_tracked` and, over a SET of teachers,
`censusAngleMapJ_foldFree_on_of_tracked` — the shape the schema actually
consumes.  The continuation certificates discharge the SECOND conjunct at the
teachers of a grid.  The first is the residual, and stating it over `p ∈ S`
rather than over a finite list is what makes visible that a certificate
visiting finitely many teachers does not meet it.  Those three are axiom-clean
`[propext, Classical.choice, Quot.sound]` and are BOOKKEEPING: they prove
nothing about the census map, and `censusZerosContainedIn_self` is committed
beside them to record that the factorisation is free at `L = K` and gains
nothing there.

Nor is a step a statement about the parameters between its endpoints: it is a
statement about its endpoint.  The parametrized route is the one that would
speak for a whole parameter interval, and it is capped at width `0.008` by the
solution's drift.

What remains, as an independent obligation, is the degenerate-zero residual: a
degenerate zero satisfies three equations in four unknowns (`map = 0` twice,
`det = 0` once), so it is a CURVE and costs `h^{-1}` where the dead covering
sweep's surface cost `h^{-2}`.  That is an exponent cheaper and it is still not
free, and the measured boundary of the determinant test — reachability ending
at `D* ≈ 0.40` against certified-domain floors of `0.02` and `0.005`, with
eight certified families in the gap across four faces — is unchanged by
anything in this section.

## The enclosure-form inventory

Two defects this campaign hit were the same shape: a repair applied to one
consumer of the enclosure layer and not to its siblings.  The Jacobian fix was
needed by three tests at once (mean-value, monotonicity, Krawczyk) and only the
first was checked; the mean-value row form was applied to the determinant test
and not to the angle-map test beside it.  So the consumers are listed, with the
form each uses and why.

| consumer | form it uses | status |
|---|---|---|
| `EXCL_PLAIN` (search + replay) | Cramer `sep_res` | measured EQUIVALENT to the row form — 0 boxes excluded by only one, in 1500 samples — so nothing to gain |
| centered test (search) | **row form first**, Cramer as fallback | the row form strictly DOMINATES: 317 both, 21 row-only, **0 Cramer-only**, and it is cheaper (4.4 vs 4.6 ms) |
| centered test (replay) | BOTH, either accepted | must keep the Cramer form: certificates emitted before the row form existed claimed leaves under it |
| Krawczyk (`_krawczyk_sD`, inflation, dedup) | Cramer Jacobian | conversion MEASURED and declined: 6 of 400 surviving boxes, and only 1.33x on the parametrized threshold |
| fold-freeness MAP test | plain row form, then **mean-value row form** | the mean-value form was missing here; at side `3e-2` it drops the survivor fraction from 32 % to 8 % |
| fold-freeness SCHUR test | mean-value determinant | with `β` swept as an interval, the literal form leaves 8 of the 55 families unresolved below `1e−8` and the mean-value form 6 — better, but not the "resolves all 45" this row used to claim, which was measured at THIN `β` on the 45 families known then (`row7a_band.py`; see "Two corrections to numbers this file published") |
| continuation PREDICTOR (`cont_cert.newton_correct`) | division-free row form + its Jacobian | a guess, never a verdict — but it is where the row form's SPURIOUS zero manifold `D = 0` bites, hence the gap floor in the corrector |
| continuation CERTIFY (`cont_cert.certify`) | Cramer Jacobian, through `census_cert._krawczyk_sD` | the same Krawczyk the enumeration uses, and the same knowing compromise as the row above it |
| continuation FOLD-FREENESS (`cont_cert.jac_det_theta`) | division-free row-form derivative | `det_theta = dE0/dD·dE1/ds − dE0/ds·dE1/dD` is `−separatedAngleJacDetJ`; gated by `cont_gate.py` at `3.4e-88` |

The one row that is still a knowing compromise is Krawczyk, and it is
compromise by measurement rather than by omission.

## The filter inventory

Every threshold that can decide a verdict, what relaxing it costs, and which
claims rest on it.  A "none found" here is evidence only within these.

| filter | value | decides | cost to relax | claims resting on it |
|---|---|---|---|---|
| `delta` gap floor | 0.02 / 0.005 | which families are inside the certified domain | domain grows into the collar where enclosures degrade; **no value covers `beta -> pi`** | EVERY per-witness completeness claim, bounded by `beta < pi - 2*delta` |
| `minw` | 2e-5 search, 1e-8 sweep | when a box is declared UNDECIDED | lower = more subdivision; raising it CANNOT hide anything, since undecided leaves FAIL the replay | none silently — this filter is self-announcing |
| `lat_tol` | 1e-5 | the exact-fit / separated SPLIT (not the class count) | a parametrized box wider than this in `beta` misclassifies the fit | `N = 2*(families + exact fit)` bookkeeping only |
| `fitr` exact-fit tube | 0.05 | which boxes the fold-freeness sweep skips as row 7b | shrinking it makes the sweep HARDER (the determinant is identically zero on the arc); enlarging shifts coverage from 7a to 7b | the fold-freeness sweep's scope |
| `dmax_pad` | 0.05 | how far `D` runs past `pi` | over-counts rather than misses — the safe direction | none (safe) |
| `seam` | 0.137 | where the `s`-period starts | chosen so the exact fit is interior; any generic offset works | none (checked by the certificate's own root-cell regeneration) |
| prefilter `EPS` | 1e-12-ish | order-only exclusions | returns 0 when uncertain, so it can only UNDER-kill | none (conservative by construction) |
| working precision | 200 bits | all enclosures | measured IRRELEVANT at the collar: 200 vs 800 bits identical to four digits | none |
| continuation `trust` | 0.08 | the damped Newton predictor's step cap | larger = the predictor jumps sheets at the degenerate column (measured: an undamped step of length `3.2` at `beta = 0.15`); smaller = more parameter bisections and more wall clock | none silently — a lost sheet is a REFUSED step, and a wrong box would fail the disjointness check |
| continuation `delta` | 0.02 | the corrector's gap floor | the row form vanishes IDENTICALLY on `D = 0`, so removing it lets the predictor fall onto a spurious zero manifold | the reading of the tracked sheets as separated families |
| continuation `maxjump` | 0.35 | when a certified step counts as having left its sheet | larger = a jump is only caught by the end-to-end disjointness check; smaller = more bisections | the TRACKING, never the per-box claim: the refused box is still a valid enclosure of a valid zero |
| continuation `r0` / `rmin` | 1e-3 / 1e-11 | the Krawczyk radius ladder | raising `rmin` CANNOT hide anything — an uncertified step is `UNDECIDED_MINRAD` and fails the replay | none — self-announcing |
| continuation `min_split` | 64 | how far a parameter interval may be halved | below it the leaf is undecided and the replay fails | none — self-announcing |

The one filter that has ever changed an answer is `delta`, and it is now
quantified rather than assumed.

The tree is walked in PREORDER, root cells in index order, child `0` before
child `1`.  At every node exactly one bit comes first: **`1` = internal**, its
two children following immediately in the stream; **`0` = leaf**, followed by
three bits of verdict:

| code | verdict | what the leaf claims |
|---|---|---|
| 0 | `PREFILTER` | the sign law excludes the box (soundness: the landed straddle / same-side theorems) |
| 1 | `EXCL_PLAIN` | the plain outward-rounded enclosure of `G0` or `G1` misses zero |
| 2 | `EXCL_CENTERED` | a mean-value enclosure `G(X) ⊆ G(x_c) + J(X)·(X − x_c)` misses zero, in EITHER admissible form (see below) |
| 3 | `KRAWCZYK_EMPTY` | `K(X) ∩ X = ∅` — no solution |
| 4 | `KRAWCZYK_UNIQUE` | `K(X) ⊆ int(X)` — EXACTLY one solution, enclosed |
| 5 | `KRAWCZYK_INFLATED` | the same on an inflated `Z' ⊇ Z`, so the leaf holds at most that one |
| 6 | `UNDECIDED_MINWIDTH` | no verdict, minimum width reached |
| 7 | `UNDECIDED_OTHER` | no verdict, any other reason |

Codes 6 and 7 are not distinguished further because the checker treats any
undecided leaf identically: it invalidates the completeness claim and the
replay exits non-zero.  Bits are packed MSB-first into bytes, the last byte
zero-padded, and base64-encoded as `bits_b64`; `n_nodes` gives the number of
meaningful bits so the padding can never be misread as extra tree.

**Why a checked tree is a completeness proof.**  Coverage is structural: a
node is a leaf or the union of its two children, so the leaves tile the root
cells and the root cells tile the domain.  If every leaf's verdict is valid —
provably no solution, or exactly one and enclosed — then the domain holds
exactly the enclosed families and nothing else.

**The two admissible mean-value forms, and why the second one exists.**  The
residual pair can be written two ways.  The CRAMER form `G_i` solves the two
radial rows for the student weights, which divides by
`massDet(D) = π² − phiCosJ(D)²`.  The ROW form `E_i = massDet(D)·G_i`
(`census_cert.sep_res_rowform`) never divides; the two have the same zero set
and the same sign, since `massDet > 0` off coincidence, so either may exclude
a box.

The difference is not in the residual but in its DERIVATIVE, which the
mean-value form leans on.  Each Cramer Jacobian entry divides by `massDet`
three times and then subtracts two same-sized quotients — a dependency
amplifier.  Measured at the F1 witness on 200 plateau-scale boxes that survive
every other test:

| enclosure | width relative to the true (sampled) range | |
|---|---|---|
| | median | p90 |
| Cramer `dG₀/dD` | 5041× | 294 538× |
| division-free `dE₀/dD` (`sep_res_rowform_jac`) | **27.9×** | **127.9×** |
| Cramer plain `G₀` | 512× | 4 309× |
| division-free CENTERED `E₀` | **5.8×** | **38.7×** |

Everything that consumes the Jacobian inherited the first row: the mean-value
form lost its second-order advantage, and the monotonicity test — the range in
a variable whose partial is sign-definite is attained at the endpoints — could
almost never fire (6 of 308 boxes; 10 of 400 even with the tightened
derivative, because the derivative genuinely turns there, which is why
monotonicity is the wrong consumer and the mean-value form is the right one).

**The same treatment does NOT help Krawczyk, and that is worth recording.**
Krawczyk consumes the same Cramer Jacobian, so it was the obvious next
candidate — the shared-input lesson applied a second time.  Measured: a
Krawczyk operator built on the division-free derivative resolves **6 of 400**
boxes that survive every test including the new centered form (`empty` 6,
`unique` 0), against 164 for the centered form.  The reason is structural
rather than numerical: Krawczyk's strength is contraction NEAR a solution,
while the boxes that survive are ones whose residual enclosure still straddles
zero away from any solution, and those are exclusion work, not contraction
work.  So the conversion is not made: it would widen the verdict semantics for
a 1.5 % gain.

**Measured payoff.**  Of 400 boxes at the F1 witness that survive the
prefilter, the plain enclosure, the Cramer centered form AND Krawczyk, the
division-free centered form excludes **164** — 41 % of the plateau, in one
pass.  On the same witness and the same gap floor `delta = 0.005`, against the
recorded run: live area `3.775 → 2.297` at step 5 001 and `2.908 → 1.506` at
step 10 001, with the heap 4.6× smaller at equal area.  And the ceiling is
total: of 200 plateau-scale surviving boxes, **200** have a sampled range of
`G₀` or `G₁` that misses zero, so every one of them is excludable in principle
and none of the plateau is a genuine near-solution set.

Both forms are gated before use by `validate_rowform.py`, which is a
FALSIFIER, not evidence: the identity `E_i = massDet·G_i` to `4.2e-46`, the
derivative against the audited Cramer Jacobian through
`dE/ds = massDet·dG/ds` and `dE_i/dD = (dmassDet/dD)·G_i + massDet·dG_i/dD` to
`9.6e-45`, and against central differences to `5.0e-13`.  An invalid
derivative enclosure is the failure this certificate has already paid for
once, so it is gated rather than trusted.  The replay checker accepts an
`EXCL_CENTERED` leaf if EITHER form signs a component, which keeps every
certificate emitted before the row form existed replayable.

**The boundary margin the replay also reports.**  The two solution verdicts
(`KRAWCZYK_UNIQUE`, `KRAWCZYK_INFLATED`) carry an enclosure of the solution
they assert, and the checker re-derives it.  Its `boundary` line prints the
smallest distance from any such enclosure to the four faces of the search
rectangle `s ∈ [seam, seam + 2π]`, `D ∈ [delta, π + dmax_pad]` — computed
from the OUTER endpoints of the enclosure, hence a rigorous lower bound on
the distance from the true solution to the boundary — and the replay exits
non-zero if it is not positive.  A positive margin is precisely the four
strict inequalities that the Lean interiority lemma turns into
`x ∈ interior K`, which is why row 6 of the chain below is discharged by
this line and by nothing else in the file.  The count printed there is of
solution-bearing LEAVES, before the fold and `2π` duplicates are merged; it
is not the family count.

Packing is a pure re-encoding, verified by round trip on every run
(`ivcert_pack.py`): measured at the F5 representative, 26524 leaves go from
1032.7 KB of path-keyed JSON to 22.1 KB of bitstring, 47x smaller, and the
packed file replays to the identical verdict.  **`ivcert_CERT5.json` is kept
in plain JSON as the worked example** — the one presented in full in the
webapp appendix as the exemplar of how every certificate works; the others
ship packed.

```
python3 ivcert_check.py ivcert_CERT5.json        # the worked example
python3 ivcert_pack.py  ivcert_<tag>.json        # -> ivcert_<tag>.bits.json
python3 ivcert_check.py ivcert_<tag>.bits.json   # replay, exit 0 = valid+complete
```

## The parametrized route: sound, and measured too small for a face

The enumeration has a mode that takes a BOX of teachers instead of one, and
its promise is large — a complete decomposition over a teacher box would
certify the family count constant across that box directly, without the Lean
schema.  The argument is sound, and it is worth stating in the explicit form,
because parametrization also introduces obligations the point certificates do
not have.

**Why a box of teachers is certified at once.**  Every quantity the search
evaluates is an outward-rounded interval computed from the parameter box `P`
together with the student box `Z`, so it encloses the value for every teacher
`p` in `P` at once.  Each verdict then specialises to a fixed `p`:

* *exclusion* — if the enclosure of `G₀` (or `G₁`, or either mean-value form)
  over `Z` with parameters ranging over `P` misses zero, then for each fixed
  `p` the value `G₀(p, z)` lies in that same enclosure for every `z ∈ Z`, so
  `G₀(p, ·)` has no zero on `Z`;
* *uniqueness* — every interval entering the Krawczyk operator at `P` encloses
  the corresponding quantity at `p`, so `K_P(Z) ⊇ K_p(Z)`; hence
  `K_P(Z) ⊆ int(Z)` forces `K_p(Z) ⊆ int(Z)`, and Krawczyk's theorem for the
  FIXED map `G(p, ·)` gives exactly one zero in `Z`, for each `p` separately;
* *emptiness* — likewise `K_P(Z) ∩ Z = ∅` forces `K_p(Z) ∩ Z = ∅`.

Coverage is structural as in the point case, so a parametrized run with no
undecided leaf would certify the count constant on `P`.

**Two obligations parametrization adds.**  (i) The order-only prefilter is not
one law but two — positive masses force each student's offsets to STRADDLE the
gap lattice, mixed masses force them to the SAME side — so the branch must be
derived from the mass interval over the whole `y`-box, and a box straddling
`y = 0` has no single law and must be refused.  Taking the branch from one
endpoint can re-verify a `PREFILTER` leaf against the wrong half of the
dichotomy: a certificate worthless while looking green.  (ii) The exact fit
sits at `θ = (β, 0)` and therefore MOVES with the box, while the dedup
identifies it by a lattice test with a fixed tolerance `1e-5`; over a wider
`β`-range the fit stops being recognised and the separated/fit split silently
goes wrong even though the class count stays right.  It must be identified by
its closed-form location intersected with the box.

**And then the measurement that decides the route.**  At the certified F4
family, parametrized Krawczyk certifies uniqueness only up to a parameter
width of about `0.006`, and `0.008` with the division-free derivative:

| parameter width | 0.006 | 0.008 | 0.0125 | 0.025 | 0.05 | 0.10 |
|---|---|---|---|---|---|---|
| Cramer Jacobian | unique | — | — | — | — | — |
| division-free | unique | unique | — | — | — | — |

A 180× tighter Jacobian buys **1.33×** in parameter width, which is what says
the binding constraint is NOT enclosure quality.  It is the solution's DRIFT:
over a parameter box the solution sweeps an arc, the enclosure must contain
the whole arc, and Krawczyk cannot contract below it.  Above the threshold a
box can never finish at all — a leaf holding a genuine solution can neither be
excluded nor certified, so it subdivides forever, which is exactly what four
scaling runs at widths `0.0125` to `0.10` did: live area falling steadily with
`sols 0` throughout.

The cost follows.  A face like F4 spans `β ∈ (0, π)` and about `0.8` in `y`;
at width `0.008` that is roughly `393 × 100 ≈ 39 000` parameter boxes, each
costing at least what a point certificate costs here (`59 233` leaves,
`1 578 s`) and in practice more.  That is upwards of a year of single-core
time for one face, and even a `0.1 × 0.1` patch needs about 156 boxes.

**So the route is sound and its scope is a small neighbourhood of a teacher,
not a face.**  The structural reading is the useful part: parameter boxes
handle the solution's drift by ENCLOSING it, which costs area in proportion to
the box width, whereas the count-constancy schema handles the same drift
ANALYTICALLY through the implicit function theorem and does not care how far
the solution moves.  That is the division of labour the chain already has —
interval certificates at points, the Lean schema across a face — and this
measurement is the reason it is the right one.  Row 7a's mass-free sweep
therefore remains the path to a face-level count.

## What is proved about this map today

One place, in three parts, so a referee does not have to assemble it from the
chain.

**Certified by interval arithmetic, and replayed** (rows 4–6a).  At the
twenty-two teachers of the table below -- every proposed label included -- the separated
family count is COMPLETE:
every point of the enumeration domain is either excluded by a certificate or
inside a Krawczyk box, so the listed families are all of them — exactly, not
"at least" — and every enclosed solution is strictly inside the search
rectangle with a positive margin. At the named replayed teachers the counts
are `3` for those assigned F2, `1` for those assigned F4, and `3` for those
assigned F6. The F5-labelled teachers have count `3` at the sampled middle
heights and count `1` at the sampled outer heights. Two F5 total-count changes
are bracketed at `beta = 3.0`, each by a pair
of complete point certificates that disagree; no wall curve is enclosed.
Each of these is a statement about the teachers
it names and about no others.

**Proved in Lean, as generic theorems and implications** (rows 1–3, 6b, 7, 8,
9, 9b).  The
enumeration solves the criticality system and nothing else; the strict sign
law of the separated scalar and the soundness of the order-only prefilter; the
enumeration domain is compact and its four strict inequalities give
interiority; fold-freeness follows from a MASS-FREE quantity; the labelled
zero count is `2·(families + exact fit)`; the count is constant on a
preconnected teacher region when the theorem's explicit boundary and
fold-freeness hypotheses hold; and the rectangle form derives finiteness
rather than assuming it.

**Not claimed: any count on a complete global F1--F7 arrangement.**  A regional
statement consumes, at every teacher of the region, fold-freeness at every
zero of the angle map.  At the
exact fit that is already certified (row 7b); away from it, it reduces to the
mass-free `generalJKernelTeacherSchurDetT ≠ 0`, whose per-face sweep does not
exist yet — though its headroom is now measured and large (smallest relative
margin `0.432` over the 55 non-exact-fit certified families, re-measured
2026-08-16, none near a fold).  The
parametrized enumeration was tested as an alternative that would need no
sweep at all, and measured too small for a proposed global cell -- sound, but limited to
parameter widths around `0.008` by the solution's drift, which no enclosure
improvement fixes. So what stands today is per-teacher completeness at
twenty-two teachers across all seven proposed labels, plus the generic Lean
propagation theorems. On the separate 68-tile rectangle, external interval
replay supplies the regional hypotheses subject to the human-checked
expression correspondence, so that one regional count is closed. Nothing
promotes it to the complete global map. F5 is the warning against that
promotion: complete point enumerations under the same proposed label have
different total family counts, even though the selected census strings at
those sampled points agree.

**What a closed face will and will not be — read this before the count.**  A
face-level count will NEVER be a hypothesis-free Lean theorem here, and it is
not one of our obligations to make it into one.  The shape of the finished
claim is fixed, and it is two halves that a referee must be able to tell apart:

* **kernel-checked half.**  `censusAngleMapJ_familyCount_eq_witness_on_censusBox`
  is a Lean theorem carrying fold-freeness as an EXPLICIT hypothesis (`hfold`,
  the Jacobian determinant nonzero at every zero over the teacher rectangle),
  alongside boundary containment.  It is machine-checked and axiom-clean, and
  it says: GIVEN those inputs, the family count at the witness is the count at
  every teacher of the rectangle.
* **replayed half.**  An interval certificate discharges `hfold` FOR THE
  CONCRETE CENSUS MAP on the named rectangle — never as a general theorem, per
  the owner's scope rule — emitted in the same bisection-tree format, replayed
  by the same checker, undecided leaves failing the replay, and presented in
  full rather than as a verdict.

The count therefore rests on a Lean implication plus a replayable computation,
and this file will say so wherever the count is stated.  It will not say the
face count is "proved in Lean": the implication is, its hypothesis is not.
That division is the point of the interval-certificate rule, not a concession
in it.

**And the continuation certificate is NOT that replayed half.**  It supplies
`separatedAngleJacDetJ ≠ 0` at the zeros it TRACKS, at the teachers of a grid.
The schema's `hfold` asks for it at every `x ∈ censusSquareJ seam delta` and
every `p` in the rectangle.  The two differ by exactly the obligation the
continuation cannot carry — that the tracked list is all of them — and that
obligation is what the per-teacher enumeration supplies at ONE teacher and
nothing supplies across a face.  See the continuation section for the full
statement of the gap.

**The remainder, with its obstruction named.** The F1 witness enumeration
closed on 2026-08-10 (the row above; 509 533 steps, zero undecided area) — the mixed sector's weak
order-only prefilter cost it a factor in steps, not a wall.  The mixed collar
is still open in the same sector.  F3 and F7 have per-teacher certificates at
`(3.08, −0.08)` and `(3.08, −0.92)`; their continuation certificates need the
gap floor `delta = 0.002`, because on those proposed-label grids the smallest family gap
scales as `c(beta)·|y|` with `c → 0` at `beta = pi`.  The columns `beta = 0` and `beta = pi` are
exact degenerate cases treated in closed form, not gaps.

## The one input the face statement is still missing, and its headroom

> **READ THE LATER TREATMENT FIRST.**  This section is the OLDEST of four
> passages in this file about the same quantity, and one of the later ones
> is a retraction.  Its box-side table is superseded and marked as such
> below.  The current statement of what the mass-free route can and cannot
> reach is "Row 7a in Lean as a REGION statement", and the number that
> decides an individual family is in "The law that says which families are
> cheap, before any of them is swept".

Row 7's remaining obligation is a single nonvanishing: at every zero of the
census angle map over a face, `generalJKernelTeacherSchurDetT ≠ 0`.  It is
MASS-FREE — its arguments are the teacher gap and the two student angles, not
the masses — which is why one certificate can serve a whole face, whose `y`
coordinate only moves the masses.  Before building that sweep, two things were
measured (`kernel_schur.py`, `kernel_schur_headroom.py`; the transcription is
gated against the tree-verified row form at `8.2e-35`, worst over the 55
families, re-run 2026-08-16).

**The exact fit is outside this route, and by a theorem rather than an
accident.**  At `θ = (β, 0)` the whole `2 × 2` angle matrix vanishes for every
teacher (`generalJAngleEq_eq_zero_at_exactFit`), so BOTH unit-teacher values
vanish, the kernel vector this function substitutes is `(0, 0)`, and a form
homogeneous of degree four in it is identically zero — measured `1e-178` and
below, which is zero at any working precision.  So the mass-free hypothesis is
UNSATISFIABLE at the one zero every teacher has, and a face certificate built
only from it would have been unusable.  The schema's other certified route
covers exactly that point: the mass-CARRYING `det T ≠ 0`, which the census
already certifies through its interval Schur type — `spurious`, a definite
form, at the exact-fit class of all 22 replayed witnesses (one such class per
witness, counted 2026-08-16).  Hence rows 7a and
7b.

**Row 7a is an EXCLUSION problem, and it has no drift wall.**  The obligation
is the contrapositive: there is no `(p, x)` in the region at which the angle
map vanishes AND the determinant vanishes.  Per box, two ONE-SIDED
alternatives discharge it — the determinant's enclosure misses zero (a test on
`(β, θ₀, θ₁)` alone, so one verdict serves EVERY `y`), or some angle-map
component's enclosure misses zero.  Nothing here encloses a moving solution or
contracts onto one, so the parameter-width threshold that limits the
parametrized route does not apply.

**But the cost is set by the ENCLOSURE, not by the margin, and that is the
finding.**  At each certified family the determinant's enclosure width grows
linearly in the box side, and the slope is family-dependent, so the test fires
only below a family-dependent box side.  Re-measured 2026-08-16 with
`kernel_schur.py` at the four F4 families, over box half-widths `1e-3` down to
`1e-6` in each of `β, θ₀, θ₁`, the width/side ratio is constant to four digits
and runs `9.5e4 × side` at `(0.30, −0.50)`, where the determinant is `−10.5`,
up to `2.7e12 × side` at `(2.00, −0.30)`, where it is `−5.2e9`; at the F4
teacher of this file's own F4 rectangle, `(0.75, −0.40)`, it is `4.1e8 × side`
against a value of `−3.6e5`.  (An earlier draft of this line reported
`9.5e7 × side` against a value of `−1541`; that pair reproduces at no F4
family and is withdrawn.)

> **SUPERSEDED TABLE — DO NOT BRIEF FROM IT**
>
> **This table is retired as of 2026-08-10.  It has cost one lane a wrong
> brief and it is left here only because it was published.**  It was measured
> at THIN `β`, on the 45 families known then, with the LITERAL enclosure, and
> every one of those three is now out of date.  Three later treatments of the
> same quantity live further down this file and one of them is a retraction; a
> reader who stops here gets the oldest of the four.  Read instead, in order:
> "The mean-value form removes the tail" (same section, below), the two
> corrections in "Row 7a in Lean as a REGION statement", and — for the number
> that actually decides a family — **"The law that says which families are
> cheap, before any of them is swept"**.
>
> | required box side | `1e-2` | `1e-3` | `1e-4` | `1e-5` | `1e-6` | `< 1e-7` |
> |---|---|---|---|---|---|---|
> | families (SUPERSEDED) | 13 | 5 | 8 | 6 | 7 | 5 |
>
> What replaced each of its claims: the "tail of five" is **six** on the 55
> families now certified, once `β` is swept as an interval rather than held
> thin; the "thirteen easy families are nearly free" costing used `L / side`,
> which counts a CURVE, while the zero set is a two-dimensional SURFACE in
> `(β, y, θ₀, θ₁)` and the measured cost is `~1.5e9` nodes per face and gap
> band; and "the tail alone would need `1e7` boxes" understates it, since the
> tail is not reachable at any box side tried and is now REFUSED by a measured
> rule rather than merely priced.

**A correction: the sweep's collar failure was a DESIGN flaw of mine, not a
degeneracy of the mathematics.**  An earlier version of this section reported
that the sweep could not discharge fold-freeness near the coincidence collar
and blamed a structural degeneracy.  Half of that was wrong and the operative
half was mine.

What is true: the determinant DOES degenerate as the gap closes.  Along `D` at
fixed `s`, the kernel vector norm falls from `8.7` at `D = 2` to `5.6e-5` at
`D = 0.02`, and `det T / ‖S‖⁴` still falls seventeen orders over that range —
the separated stratum merging into the coincident one as the two students
approach.  So the SCHUR test genuinely cannot fire at the collar.

What was wrong: that this blocks the sweep.  It does not, because the MAP test
covers the collar — once the sweep is built correctly.  The first version swept
only `(β, s, D)` and passed `y` as a FIXED interval, so the `y`-variation put a
floor under the angle-map enclosure that no `(β, s, D)` subdivision could beat.
Measured at the collar corner, at a box side of `1e-6`:

| `y`-width | 2e-3 | 2e-4 | 2e-5 | 2e-6 | 0 |
|---|---|---|---|---|---|
| residual relative width | 127 | 13.0 | **1.52** | 0.37 | 0.25 |
| MAP fires | no | no | **yes** | yes | yes |

The obligation is four-dimensional and so is the sweep.  With `y` subdivided,
the same cell that produced ZERO kills in 4 067 nodes gives **8 875 MAP kills
in 17 796 nodes**, at `11.2` ms per node instead of `29.5`, with the stack down
to 47 — nearly exhausted.  The determinant test staying `y`-free remains a
genuine advantage: one SCHUR verdict still serves every `y` in a box.

Two false leads were checked and cleared on the way, and both are worth
recording because they look plausible: raising the working precision from 200
to 800 bits changes nothing (identical to four digits), and NORMALISING the
vanishing factors cannot help at all — an interval sign test is
scale-equivariant, so dividing by a positive constant leaves `width/|value|`
exactly unchanged.  The reason the division-free row form and the mean-value
form did work is that each removed a dependency-amplifying OPERATION, not a
scale.

**The mean-value form removes the tail, and the sweep becomes costable.**
Built as `kernel_schur_mv.py`: the determinant is a composition eight deep
(`num` → kernel vector → `V, L, N₀, N₁` → `T₀₀, T₁₁, T₀₁` → `det`), so its
derivative is taken by FORWARD-MODE AUTOMATIC DIFFERENTIATION over intervals —
each quantity carried as `(value, ∂/∂s, ∂/∂D)` — rather than by hand, because
an eight-deep chain is where a sign error hides.  Gated three ways before use
(`kernel_schur_mv_gate.py`; re-run 2026-08-16, all four gates green): the AD
value reproduces the literal determinant to `3.3e-44` over 200 samples, the AD
derivative matches central differences to `5.6e-14`, and the
mean-value interval contains the true value at 3 000 samples with none outside.
The gate earned itself immediately: the first draft had `d|sin t|/dt = +cos p`
where the tree's branch coordinate `p = π(1+2⌊n/2⌋) − t` has `dp/dt = −1` on
EVERY branch, so the correct seed is `−cos p`.  The atoms' own `dphi` and `dsA`
already carry that convention, which is why they checked out while the
hand-seeded `|sin|′` did not.

Required box side across the same 45 families:

| | resolved | largest | median | smallest | unresolved |
|---|---|---|---|---|---|
| literal enclosure | 40 / 45 | `5e-2` | `2e-4` | `2e-7` | **5** |
| mean-value form | **45 / 45** | `1e-1` | **`2e-3`** | `5e-7` | **0** |

The five-family tail that was the whole cost is gone — those families gain
between 20× and 1000× — and the median improves tenfold.

> **CORRECTED 2026-08-10: that table was measured at THIN `β`, and a face sweep
> does not sweep `β` thinly.**  Re-measured with `β` carried as an interval of
> the same width as the box side, over the 55 non-exact-fit families now in the
> replayed witnesses (`row7a_band.py`), the mean-value form resolves 36 of 55 at
> a side of `1e−4` or coarser and leaves **six** unresolved below `1e−8` — the
> tail is smaller than the literal form's eight, not empty.  The per-face
> breakdown is in "Row 7a in Lean as a REGION statement" below; the short form
> is that F1 and F4 are the faces all of whose families the determinant test can
> serve at a practical side, and F3 and F7 are the faces where it serves none.

Costing the THIN-`β` sweep of the table above from `Σ L / side` over those
45 families, a face of `β`-range `π` needs about
`9e6` boxes, dominated by the single worst family at `5e-7`, against an
unbounded figure before (the tail was unresolved at `1e-7`, so it could not be
costed at all).  Twenty-eight of those 45 need only `1e-3` or coarser.  Both
figures are properties of the 45-family thin-`β` measurement, not of the
55-family swept-`β` one that corrects it.

Note what this does NOT follow from: every one of the 55 non-exact-fit
families has a relative margin of at least `0.432`, so the POINT margin — large
everywhere — does not predict the box cost at all.  What sets the cost is the
enclosure's dependency slope.  The repair class is therefore the one that already worked
twice here: a mean-value form for the determinant built on its derivative.
The obvious algebraic cancellation was checked first and is not the culprit —
`detT = T₀₀T₁₁ − T₀₁²` has an exactly cancelling `A²` term whose removal is
verified to `2.2e-44` and buys a median enclosure factor of `0.97`, i.e.
nothing (re-measured 2026-08-16 at all 55 families, box half-width `1e-4`;
the ratio runs `0.30` to `1.44` across them).

**Everywhere else the headroom is large.**  Re-measured 2026-08-16 by
`python3 kernel_schur_headroom.py` over the **55** non-exact-fit certified
families read from the 22 replayed witnesses (`replay_C_*.log`, all 22 carrying
`CERTIFICATE VALID AND COMPLETE`; the 22 exact-fit classes, one per witness,
are excluded for the reason above), the relative margin
`|det T| / max(|T₀₀T₁₁|, |T₀₁²|)` — how far the determinant stays from the
cancellation of the two terms that form it — runs:

| | smallest | median | largest |
|---|---|---|---|
| relative margin | **0.432** (F5 at `(3.00, −0.42)` and its `y`-mirror at `(3.00, −0.58)`, `D = 1.9373`) | **1.000** | **1.813** (the two F5 families at `(2.40, −0.50)`, `D = 2.4166`) |

Twenty-nine of the 55 lie between `0.61` and `1.17`; four are below `0.5`;
and **zero** are below `1e-6` — the closest approach to that floor clears it by
a factor of `4.3e5`.  Nor does the margin track the gap: the two smallest-gap
families (F3 and F7 at `β = 3.08`, `D = 0.0315`) sit at `0.610`, well above the
minimum, which is at `D = 1.94`. None of these 55 named witnesses is near a
determinant zero. This supports attempting a regional sweep; it is not a
fold-freeness statement for any complete proposed face.

> **CORRECTED 2026-08-16.**  The table above previously quantified over "the 42
> non-exact-fit certified families" and read `0.205` (F6 at `(1.5, −0.96)`,
> `D = 0.358`) / typical `0.6 – 0.9` / `1.007`.  Every one of those four figures
> is retired.  The set is not 42 — the replayed witnesses hold 55 non-exact-fit
> families, the same count this file already corrected `45 → 55` elsewhere —
> and none of the three extremes reproduces under the instrument the section
> itself names: the F6 family at `(1.5, −0.96)`, `D = 0.3581`, now measures
> `0.635` and is only the eighteenth smallest of the 55, and just 14 of the 55
> fall in the old `0.6 – 0.9` band.  What survived the re-run unchanged is the
> load-bearing part: **zero** families below `1e-6`, and nothing near a fold.

## Row 7a in Lean as a REGION statement, and the three holes the mass-free route has

Row 7a's Lean half was a statement about ONE zero.
`FoldFreeInputJ.censusAngleMapJ_foldFree_of_kernelTeacherSchurDetT` takes a zero
of the census angle map, four side conditions, and
`generalJKernelTeacherSchurDetT ≠ 0`, and returns the Jacobian nondegeneracy the
counting schema wants.  What the schema actually consumes is `CensusFoldFreeOn
S K` — over a whole teacher set and a whole student region.
`CountConstancyJ/SchurExclusionJ.lean` is that lift, and the lift is where the
side conditions are paid for rather than assumed.  Twenty-one declarations —
five predicates and sixteen theorems — every theorem axiom-clean at
`[propext, Classical.choice, Quot.sound]`.

**The obligation stated over a region containing the exact fit is FALSE — and
that is now a theorem rather than a paragraph.**
`generalJKernelTeacherSchurDetT_eq_zero_at_exactFitSwap` proves
`generalJKernelTeacherSchurDetT β ![β, 0] = 0` from two landed facts: both
unit-teacher values of the first eliminated angle equation vanish at `θ = (β,0)`
(`generalJAngleEq_eq_zero_at_exactFit_swap`), so the substituted kernel vector
is `(0,0)`, and the cleared Schur determinant is homogeneous of degree four in
that vector (`generalJCramerSchurDetT_smul`).  Since the exact fit is a zero of
the census angle map at EVERY teacher and lies in the enumeration domain,
`not_censusSchurExclusionOn_of_exactFit_mem` refutes the naive exclusion
outright.  A certificate campaign aimed at it would have been chasing a target
that does not exist; `CensusSchurExclusionOffFitOn`, with that one
teacher-dependent point removed, is the one an interval sweep can attempt.

**Two of the four side conditions are free on the enumeration domain, and the
gap ceiling is SHARP.**  `censusStripJ seam delta dmax` pins `θ₀ − θ₁` into
`[delta, dmax]`, so `0 < delta` excludes coincidence and `dmax < π` excludes
opposition — `not_coincident_of_mem_censusStripJ`,
`not_opposite_of_mem_censusStripJ`, derived, not assumed.  The ceiling cannot be
relaxed: `exists_opposite_mem_censusStripJ` exhibits an exactly antipodal pair
in the domain as soon as `dmax ≥ π`, and this is not an academic corner —

> **five of the 55 certified non-exact-fit families sit at `D = π` EXACTLY**, all
> of them on the equal-teacher-mass line `y = −1/2`: `F4 (0.30, −0.50)` and
> `F5 (2.40, −0.50)`, `(2.60, −0.50)`, `(2.80, −0.50)`, `(3.00, −0.50)`, each
> with `D = 3.141592654` in its replay log.

Those five are outside `GeneralJSeparatedPairCritical` by that predicate's own
definition — non-opposition is one of its FIELDS, not a lemma about it — so the
mass-free route cannot reach them for a definitional reason rather than a
numerical one.  The labelled square is worse: `exists_opposite_mem_censusSquareJ`
shows it contains antipodal pairs for any gap floor `delta ≤ π`, so the strip is
the domain this route can serve and the square is not.

**What is left is exactly two interval-testable predicates**, plus one input
that is often free — and one of the two now has a cheaper, differently-shaped
alternative.  `censusAngleMapJ_foldFreeOn_censusStripJ` delivers
`CensusFoldFreeOn S (censusStripJ seam delta dmax)` from
`CensusSchurExclusionOffFitOn`, `CensusTorqueAliveOffFitOn`, and
`CensusFoldFreeAtExactFitOn` — and the third costs nothing when the gap window
starts above every teacher gap in `S`
(`censusFoldFreeAtExactFitOn_censusStripJ_of_lt_delta`), because then the exact
fit is not in the region at all.  The torque conditions are NOT free in general:
the tree's only region-shaped torque tool
(`generalTorqueJ_zero_mem_open_cells`) needs both teacher masses positive, which
the census carries only in the `y < 0` sector — and the tree had no positivity
lemma for `censusMassesJ` at all, only `censusMassesJ_ne_zero`, so that tool
could not be pointed at the census map.  `censusMassesJ_pos_of_mem_Ioo` supplies
it, and `censusTorqueAliveOffFitOn_of_studentsOffCells` then trades the torque
enclosures for `CensusStudentsOffTorqueCellsOn` — each student sitting, modulo a
full turn, outside the two open cells `(0, β)` and `(π, π + β)` that the teacher
atoms cut.  A geometric obligation instead of an interval one, stated and not
discharged, and the shape a sign-law or arc argument can reach.

### The exclusion certificates, and exactly what they cover

`row7a_excl.py` emits and `row7a_check.py` replays a certificate whose two
verdicts discharge those two predicates, in the same replayable bisection-tree
format as the rest of this directory: root box regenerated from the spec alone,
deterministic split (widest side, child 0 the lower half), preorder bits, `1`
internal and `0` leaf plus two verdict bits, an undecided leaf FAILING the
replay.  The checker also refuses any spec that does not meet the Lean
statement's own conditions (`0 < delta`, `dmax < π`, `b1 < delta`).

The two verdicts, and why each is one-sided:

* **MAP** — some component of the angle map has an enclosure missing zero over
  the whole box, so the map cannot vanish anywhere in it.  The division-free row
  form is used, and `E_i = massDet(D)·G_i` (identity re-gated this round at
  `4.2e−46` by `validate_rowform.py`), so `E_i ≠ 0` gives `G_i ≠ 0`, which is the
  direction needed.  That the PAIR really is the census angle map's pair is
  gated too, and the second row was not gated before: `row7a_rowgate.py`
  reassembles both rows from the transcribed `separatedNumJ` at the arguments
  the tree's unit-teacher lemmas use and matches `sep_res_rowform` to `4.6e−46`
  and `3.6e−46` respectively, up to a global sign that cannot change whether an
  enclosure misses zero.
* **SAFE** — the mean-value enclosure of `generalJKernelTeacherSchurDetT` misses
  zero AND both teacher-torque enclosures miss zero, so wherever the map does
  vanish in that box, the mass-free route applies and concludes.

Coverage is structural: the two children of a node are the two halves about the
midpoint of its widest side, and the `s`-range is built from the UPPER bound of
`2π`, so the root over-covers the strip rather than leaving a sliver.

Six historical artifacts are committed, but only the map-only smoke artifact
still replays after the 2026-08-10 correction to the odd-half-branch formula
for `|sin|`.  The former evaluator used `sin(p)` instead of `e*sin(p)` on odd
branches.  That does not affect a `MAP` leaf, but it invalidates every stored
`SAFE` leaf.  The checker now rejects all five artifacts containing such
leaves; they are retained as regression fixtures, not certificates.

| certificate | teacher rectangle | gap band | leaves | MAP | SAFE | family inside | replay |
|---|---|---|---|---|---|---|---|
| `row7a_smoke` | `[0.7500,0.7510] × [−0.4000,−0.3990]` | `[2.00, 2.10]` | 26 | 26 | 0 | none | VALID AND COMPLETE |
| `row7a_safesmoke` | `[0.7500,0.7501] × [−0.0201,−0.0200]` | `[1.40, 1.50]` | 399 | 257 | 142 | F2's `D = 1.4530` | **INVALID: 142 mismatches** |
| `row7a_scale0.0002` | `[0.7500,0.7502] × [−0.0201,−0.0199]` | `[1.40, 1.50]` | 416 | 264 | 152 | the same | **INVALID: 152 mismatches** |
| `row7a_scale0.0005` | `[0.7500,0.7505] × [−0.0201,−0.0196]` | `[1.40, 1.50]` | 555 | 326 | 229 | the same | **INVALID: 229 mismatches** |
| `row7a_scale0.001` | `[0.7500,0.7510] × [−0.0201,−0.0191]` | `[1.40, 1.50]` | 1 455 | 895 | 560 | the same | **INVALID: 560 mismatches** |
| `row7a_scale0.002` | `[0.7500,0.7520] × [−0.0201,−0.0181]` | `[1.40, 1.50]` | 6 117 | 3 734 | 2 383 | the same | **INVALID after evaluator correction** |

What the one valid row licenses, exactly: for every teacher in that rectangle and every
student pair in `censusStripJ 0.137 delta dmax` with that band's `delta` and
`dmax`, the census angle map has no zero at which the mass-free type-boundary
function vanishes, and both teacher torques are alive at every zero — i.e.
`CensusSchurExclusionOffFitOn` and `CensusTorqueAliveOffFitOn` hold there, and
since every band starts above the teacher gap, so does
`CensusFoldFreeAtExactFitOn`.  Feeding those to
`censusAngleMapJ_foldFreeOn_censusStripJ` gives `CensusFoldFreeOn` on that
region — the schema's `hfold`, over a region, from a replayed computation and a
machine-checked implication, with the two halves separate exactly as this file's
division-of-labour paragraph requires.

The certificate proves slightly MORE than the predicate needs, and knowingly:
the torque enclosures are required to miss zero on the whole SAFE box, whereas
the predicate asks only for the zeros of the map.  That is the conservative
direction, and it is what makes the torque test cost real nodes — a box that
straddles a torque zero cannot be SAFE and has to be subdivided until the MAP
test fires.  Refining it to "no map zero inside the torque-dead part of the box"
is the obvious cheapening and is not implemented.

What no row licenses: anything about the rest of the face, anything about gaps
below the band, and anything about the labelled square, whose antipodal pairs
this route cannot serve at all.

Scope, in the terms the ledger uses: each row is a statement about the CONCRETE
census map on that teacher rectangle and that gap band, and about nothing else,
and never the basis of a general theorem.

### Two corrections to numbers this file published

**(a) The required box side was measured at THIN `β`, and a face sweep does not
sweep `β` thinly.**  Re-measured over the 55 non-exact-fit families now in the
replayed witnesses (`row7a_band.py`, `row7a_band.log`), with `β` carried as an
interval of the same width as the box side:

| form | resolved at side ≥ `1e−4` | unresolved below `1e−8` |
|---|---|---|
| literal, `β` thin | 29 / 55 | 8 |
| mean-value, `β` thin | 44 / 55 | 1 |
| literal, `β` swept | 28 / 55 | 8 |
| **mean-value, `β` swept** | **36 / 55** | **6** |

So the mean-value form's headline — the tail is gone — holds at thin `β` and
does NOT survive `β` being swept: six families stay unresolved below `1e−8`.
Per face — and the scope matters: these are the families certified at that
face's REPLAYED WITNESSES, not at every teacher of the face — resolvable at a
side of `1e−4` or coarser with `β` swept:

| face | F1 | F2 | F3 | F4 | F5 | F6 | F7 |
|---|---|---|---|---|---|---|---|
| families | 4 | 12 | 3 | 4 | 20 | 9 | 3 |
| resolvable at ≥ `1e−4` | **4** | 8 | **0** | **4** | 17 | 3 | **0** |
| unresolved below `1e−8` | 0 | 0 | **3** | 0 | 0 | 1 | **2** |

F1 and F4 are the faces every one of whose witness families the determinant test
can serve at a practical box side; F3 and F7 — the two faces at `β ≈ 3.08`,
where every family sits at `D ≤ 0.10` — are the faces it can serve at NO box
side tried, down to `1e−8`.  This ranks the faces by how hard the mass-free
route finds them; it does not close any of them, because the covering cost is
what decides that and it is measured in the next section.

**(b) and (c) The small-gap wall is DEPENDENCY, not a vanishing quantity.**
Both corrections are recorded where the claim sits — in "Why the fold-freeness
sweep does not scale to a face", above — because that is where a reader meets
them.  In one line: the linearised distance to a zero of the determinant at the
smallest-gap families is `1.0e−3`, not `1e−25`, confirmed by two computations
that share no code (central differences of the literal determinant, and the
determinant's sampled true range over a box); the determinant varies by 10 % of
its value across a box of side `1e−4` there, while the enclosure is `3e21` times
wider than that variation.  An exact-range evaluation would sign every certified
family at a side of `1e−4`.  What is NOT changed by this is the verdict: the
covering cost is dimensional, so a better enclosure buys the small-gap band and
still does not buy a face.

### The law that says which families are cheap, before any of them is swept

The tail is not merely expensive; it is predictable, and cheaply.  Two
measurements chain into a rule that costs one interval evaluation and one
sample grid per family — under a second — and decides whether a covering
certificate can ever reach that family, without running a bisection at all.

**Step 1: the gap predicts the enclosure.**  Define the DEPENDENCY FACTOR of a
family as

> `F = (width of the determinant's enclosure over a box) / (the determinant's
> true range over that same box)`,

at a fixed box side, with the true range taken by dense sampling.  `F = 1` is an
exact-range evaluation; `F = 10³` means the interval arithmetic has thrown away
three decades to dependency.  Measured at side `1e−4` over the 55 non-exact-fit
certified families (`row7a_stagewidth.py 160 all` → `row7a_depfactor.log`):

| `D` band | `0.00–0.05` | `0.05–0.15` | `0.15–0.30` | `0.30–0.50` | `0.50–1.00` | `1.00–2.00` | `2.00–4.00` |
|---|---|---|---|---|---|---|---|
| families | 2 | 4 | 3 | 8 | 5 | 15 | 18 |
| `F` | `3e21`–`7e21` | `2e14`–`3e18` | `7e05`–`8e08` | `1e02`–`2e05` | `9e00`–`4e02` | `2e00`–`1e03` | `1e00`–`3e02` |

**`F ~ D^{−9.17}`**, a log-log fit over four decades of student gap, with `2.3`
decades of residual scatter.  The scatter is real and is why the law is stated
in two steps rather than as a formula in `D`: `F6 (0.75,−0.98)` carries
`F = 1.4e+03` at a comfortable gap of `1.4530`, so a family can be expensive
without being narrow.  `D` predicts `F` on average; `F` is what actually decides.

**Step 2: the enclosure predicts the cost.**  Against the measured required box
side of the same families (`row7a_band.log`, mean-value form, `β` swept), over
the 51 families the two tables share — four are lost to key collisions, families
of one witness with the same gap to four decimals —

> Spearman `ρ = −0.910` between `log F` and `log (required box side)`,

and the decision rule that falls out of it is **one-sided and exception-free in
the direction that matters**:

| rule | families | verdict |
|---|---|---|
| `F ≥ 1e3` at side `1e−4` | 15 | **0 of 15** are resolvable at a box side `≥ 1e−4`.  No exception at any threshold from `1e3` to `1e6` |
| `F < 1e3` | 36 | 32 of 36 are.  The four failures are `F5 (3.00,−0.60)`, `F6 (0.75,−0.98)`, `F6 (2.2,−0.94)`, `F2 (1.5,−0.04)`, all with `F` in `[2.8e2, 9.4e2]` — the top decade below the threshold — and all missing by exactly one decade, needing `1e−5` |

So `F ≥ 1e3` is a REFUSAL certificate for a family, and it never lied on this
data; `F < 1e3` is an expectation of affordability that is right 32 times in 36
and wrong only just below the threshold and only by one decade.  All six
families that no box side down to `1e−8` resolves carry `F ≥ 8.3e+08`.

**Validity domain**, since the rule is only as good as its instrument.  Box side
`1e−4`; 160-bit interval arithmetic; the mean-value form of
`generalJKernelTeacherSchurDetT` with `β` carried as an interval of the same
width; the 55 non-exact-fit families certified at the replayed witnesses, and
nothing about any other teacher; the true range taken by a `9 × 9 × 3` sample
grid, which UNDER-estimates a range (§1p) and therefore OVER-estimates `F` —
the error direction that makes a refusal conservative rather than optimistic.
The rule is about the determinant test alone; it says nothing about the MAP
test, which is what covers everything away from a family.

**What this converts.**  Before, the five- (now six-) family tail was a cost:
somewhere below `1e−7` per family, priced at `1e7` boxes and rising.  Now it is
a property with a number attached — those families sit where four decades of
gap have cost twenty-one decades of enclosure — and any successor can sort a new
face's families into reachable and unreachable in seconds, rather than by
discovering it in a sweep that does not terminate.

### The price of a face, measured on the certificate rather than argued

Holding the gap band fixed at `[1.40, 1.50]` (one certified family inside) and
growing the teacher rectangle:

| teacher width `w` (both `β` and `y`) | `1e−4` | `2e−4` | `5e−4` | `1e−3` | `2e−3` |
|---|---|---|---|---|---|
| nodes to completion | 797 | 831 | 1 109 | 2 909 | 12 233 |

Below `w ≈ 1e−3` the cost is nearly flat in `w` — a factor `3.6` over a
DECADE — because the family barely moves across the rectangle and one column of
boxes covers the whole teacher extent.  From `1e−3` to `2e−3` it is `4.2×` for a
doubling of both sides, i.e. exactly proportional to teacher AREA: the family's
drift has passed the box side the SAFE test needs, and from there the cost is
the covering of a two-dimensional zero surface.  That is what makes a face
unaffordable, and now with a measured constant rather than a dimension count:
`12 233 × (1/2e−3) × (0.5/2e−3) ≈ 1.5e9` nodes for a face of `β`-range `1` and
`y`-range `0.5`, for ONE gap band of width `0.1`, at the measured `19` ms per
node — about a year of single-core time.  That is the same verdict
this file already recorded for the covering sweep, now with a measured exponent
and a measured constant behind it rather than a dimension count.

The contrast that isolates the small-gap wall inside the sweep itself, same
teacher rectangle `[0.75,0.76] × [−0.40,−0.39]`, same code, only the gap floor
changed:

| gap band | nodes | MAP kills | SCHUR kills | outcome |
|---|---|---|---|---|
| `[1.00, 3.19]` | 55 375 | 17 792 | 9 896 | **complete** in 974 s |
| `[0.02, 3.19]` | 195 000+ | 97 478 | **0** | 3 000 s timeout (`EXIT 124`), still open |

Zero SCHUR kills in 195 000 nodes is the small-gap wall as the sweep sees
it: below `D ≈ 0.4` the determinant test never fires, so every box has to be
killed by the map test alone, and the boxes that hold a zero cannot be.

### The accounting: what row 7a closes, and what it does not

**Closed after the odd-branch correction.**  Only the map-only smoke row
replays.  The five historical SAFE rows fail by exactly the leaves that relied
on the incorrect odd-branch absolute sine and license no statement:

| certificate | region | leaves | MAP | SAFE |
|---|---|---|---|---|
| `row7a_smoke` | `[0.7500,0.7510] × [−0.4000,−0.3990]`, gap `[2.00,2.10]` | 26 | 26 | 0 |
| `row7a_safesmoke` | `[0.7500,0.7501] × [−0.0201,−0.0200]`, gap `[1.40,1.50]` | 399 | 257 | **142 invalid** |
| `row7a_scale0.0002` | `[0.7500,0.7502] × [−0.0201,−0.0199]`, gap `[1.40,1.50]` | 416 | 264 | **152 invalid** |
| `row7a_scale0.0005` | `[0.7500,0.7505] × [−0.0201,−0.0196]`, gap `[1.40,1.50]` | 555 | 326 | **229 invalid** |
| `row7a_scale0.001` | `[0.7500,0.7510] × [−0.0201,−0.0191]`, gap `[1.40,1.50]` | 1 455 | 895 | **560 invalid** |
| `row7a_scale0.002` | `[0.7500,0.7520] × [−0.0201,−0.0181]`, gap `[1.40,1.50]` | 6 117 | 3 734 | **2 383 invalid** |

Only `row7a_smoke` discharges `CensusSchurExclusionOffFitOn` and
`CensusTorqueAliveOffFitOn` on its region, and its band starts above the teacher gap so
`CensusFoldFreeAtExactFitOn` is free there; through
`censusAngleMapJ_foldFreeOn_censusStripJ` that is `CensusFoldFreeOn` — the
schema's `hfold` — on that region.  The map test alone cleared every leaf, so that band provably
holds no zero of the census angle map at all, and its exclusion is true for the
stronger reason.

**NOT closed: no face, and the two attempts are priced rather than pending.**
Both runs at a `0.01 × 0.01` teacher rectangle over the full band `[1.00, 3.10]`
were still running when they were stopped, with no undecided leaves but no
completion:

| run | region | nodes reached | MAP | SAFE | state |
|---|---|---|---|---|---|
| `row7a_F2band` | `[0.75,0.76] × [−0.03,−0.02]` | 220 000 | 90 442 | 19 552 | stopped at 3 506 s, `0.612` of the volume decided |
| `row7a_F4band` | `[0.75,0.76] × [−0.40,−0.39]` | 205 000 | 69 599 | 32 884 | stopped at 3 364 s, `0.671` decided |

They are not restarted, because the scaling measurement below already says what
they would cost: a `0.01 × 0.01` rectangle is `25×` the area of the largest
certificate above and its band is `21×` wider, so both are `1e6`-node runs, and
a face is `1.5e9`.  Neither is a certificate and neither is claimed as one.

**Per face, what the determinant test could serve if cost were free.**  Over the
families certified at each face's replayed witnesses — 51 of the 55 join the
dependency table, four lost to key collisions — with the refusal rule of the
previous section applied:

| face | families | reachable at side ≥ `1e−4` | REFUSED by `F ≥ 1e3` | unresolved to `1e−8` |
|---|---|---|---|---|
| F1 | 4 | **4** | 0 | 0 |
| F2 | 12 | 8 | 3 | 0 |
| F3 | 3 | **0** | **3** | 3 |
| F4 | 4 | **4** | 0 | 0 |
| F5 | 16 | 13 | 2 | 0 |
| F6 | 9 | 3 | 4 | 1 |
| F7 | 3 | **0** | **3** | 2 |
| all | 51 | 32 | 15 | 6 |

Read this correctly, because it is the answer to "which faces do the cheap
families close" and the answer is *none, and two of them cannot be closed by
this route at all*:

* **F1 and F4** are the only faces every one of whose witness families the
  determinant test can serve at a practical box side.  They are still not
  closed — closing them is the `1.5e9`-node covering, not the per-family
  reachability — but they are the faces where the route has no in-principle
  obstruction and only a cost one.
* **F3 and F7** are refused outright: every family at their witnesses carries
  `F ≥ 2e14`, and five of the six families that no box side down to `1e−8`
  resolves are theirs.  Both sit at `β ≈ 3.08` where every family has
  `D ≤ 0.10`.  No enclosure improvement of the kind in place moves them; a
  stage-wise centered arithmetic is the only candidate and is not built.
* **F2, F5, F6** are mixed: a majority reachable, a minority refused, so a
  certificate over such a face would stall on its own worst family however
  cheap the rest are.  A covering does not get partial credit.
* Separately and independently of cost, **F4 and F5 also hold the five
  antipodal families** at `D = π` exactly, which the mass-free route cannot
  reach at any precision because non-opposition is a FIELD of
  `GeneralJSeparatedPairCritical`.  So F4's "no obstruction" above is about the
  determinant test only; the route as a whole still misses that family.

**The honest one-line summary.** Row 7a's Lean half and conditional region lift
are done. After correcting the odd-branch absolute-sine evaluator, its interval
half is discharged only by the 26-leaf map-only `row7a_smoke`; the five
SAFE-bearing artifacts are invalid regression fixtures. No proposed-label
region has a valid row-7a cover. The measured cost remains about `1.5e9`
boxes per proposed face and gap band, which is not a gap this route closes.

## The completeness chain, row by row

Owner artifact standard (2026-08-08): inside the completeness certificate the
only evidence classes are **Lean**, **iv certificates**, and this README's
glue.  Python/mpmath scripts and float logs are prospecting — they find the
numbers, they do not stand in the chain.  Each row below names its claim and
points at exactly one Lean declaration or one iv file; a referee replays the
chain from the Lean checker, `ivcert_check.py`, and this table.

| # | claim | evidence | artifact | status |
|---|---|---|---|---|
| 1 | separated criticality ⇔ the system the enumeration solves | Lean | `SeparatedCount/CensusResidualJ.lean`: `separatedPairCritical_censusResidual_eq_zero`, `separatedPairCritical_of_censusResidual_eq_zero` | **done** |
| 2 | the strict sign law of the separated scalar | Lean | `SeparatedCount/SignLawJ.lean`: `separatedNumJ_neg_of_offset_lt_gap`, `separatedNumJ_pos_of_gap_lt_offset`, `separatedNumJ_eq_zero_iff_lattice` | **done** |
| 3 | soundness of the order-only box prefilter | Lean | `LatticeStraddleJ.lean`: `separatedPairCritical_student0/1_offsets_straddle`; `MixedSameSideJ.lean`: `..._offsets_same_side` | **done** |
| 4 | per-witness completeness: the domain holds exactly the enclosed families | iv | `ivcert_<tag>.bits.json` replayed by `ivcert_check.py` | **done at the witnesses in the table below** |
| 5 | the per-witness family COUNT | iv | the same certificate | as row 4 |
| 6a | every enclosed solution is strictly inside the search rectangle, with a positive margin | iv | the `boundary` line of `replay_<tag>.log` | **done at the same witnesses**, smallest margin `5.8e-3` |
| 6b | those four strict inequalities ARE the schema's `hint`, and the schema instance on the census domain | Lean | `CountConstancyJ/CensusDomainJ.lean`: `censusAngleMapJ_zeroCount_eq_witness_strip` | **done** |
| 7a | fold-freeness from the MASS-FREE type-boundary function, away from the exact fit, at ONE zero | Lean | `CountConstancyJ/FoldFreeInputJ.lean`: `censusAngleMapJ_foldFree_of_kernelTeacherSchurDetT` | **done** |
| 7a-region | the same over a whole teacher set and the whole enumeration domain, with the lattice side conditions DERIVED from the gap window | Lean | `CountConstancyJ/SchurExclusionJ.lean`: `censusAngleMapJ_foldFreeOn_censusStripJ`, and `not_censusSchurExclusionOn_of_exactFit_mem` for why the exact fit must be removed | **done**; reduces row 7a to exactly two interval-testable predicates |
| 7a-iv | those two predicates, discharged | iv | `row7a_<tag>.json` emitted by `row7a_excl.py`, replayed by `row7a_check.py`; gated by `row7a_rowgate.py` | **one valid map-only smoke** (`row7a_smoke`, 26 MAP leaves, no SAFE leaves). `row7a_safesmoke` and four `row7a_scale*` files are invalid after the evaluator correction and license no Schur exclusion; no proposed-label region is covered |
| 7a-holes | the three places the mass-free route CANNOT reach, each named rather than left to a sweep to discover | Lean | `CountConstancyJ/SchurExclusionJ.lean`: `generalJKernelTeacherSchurDetT_eq_zero_at_exactFitSwap` (the determinant vanishes IDENTICALLY at the exact fit, so `not_censusSchurExclusionOn_of_exactFit_mem` refutes the naive obligation outright); `exists_opposite_mem_censusStripJ` / `exists_opposite_mem_censusSquareJ` (the gap ceiling `dmax < π` is SHARP, and five certified families sit at `D = π` exactly, outside `GeneralJSeparatedPairCritical` by that predicate's own non-opposition FIELD); `censusTorqueAliveOffFitOn_of_studentsOffCells` (the torque conditions REDUCED to a geometric condition, on the positivity `censusMassesJ_pos_of_mem_Ioo` the tree was missing) | **done** — hole (i) is a theorem, hole (ii) is a theorem plus five measured witnesses, hole (iii) is reduced and NOT discharged |
| 7a-gate | which families the determinant test can reach at all, decided before any sweep | **PROSPECTING** (a measurement, so not evidence in this chain — it decides where to spend, it certifies nothing) | `row7a_depfactor.log` + `row7a_band.log`, via `row7a_stagewidth.py 160 all`; the law and its validity domain are in "The law that says which families are cheap" | dependency factor `F ~ D^{−9.17}`; `F ≥ 1e3` REFUSES (0 of 15 reachable at side ≥ `1e−4`), `F < 1e3` right 32 of 36; refusals conservative because sampling under-estimates the true range |
| 7b | fold-freeness AT the exact fit, where the mass-free route cannot reach | iv | the interval Schur type of the exact-fit class in every `replay_<tag>.log` | **done at all 22 witnesses** (`spurious`, i.e. a definite form, so `det T ≠ 0`) |
| 8 | labelled-vs-unordered count, `N = 2·(families + exact fit)` | Lean | `SeparatedCount/SwapCountJ.lean`: `censusAngleMapJ_zeroCount_censusSquare_eq_two_mul_add_one` | **done** |
| 9 | count constancy across a preconnected teacher region | Lean | `CountConstancyJ/FaceCountSchemaJ.lean`: `censusAngleMapJ_zeroCount_eq_witness` | conditional schema landed; its regional `hint` and `hfold` premises are not discharged on a whole proposed label by the current row-7a artifacts |
| 9b | conditional witness-to-region family-count statement | Lean | `CountConstancyJ/SquareCountJ.lean`: `censusAngleMapJ_familyCount_eq_witness_censusSquare` | generic theorem landed; concrete F1--F7 boundary/fold hypotheses remain open, except on the separately replayed 68-tile rectangle |
| 10 | proposed-label witness lists, band agreement, collar-wall samples, F5 count band | script output | `faces_fixed.json`, `band_comparefix2.json`, `f5_count_band.json`, `witness_closure.json` | **PROSPECTING** as regional evidence; the twenty-two exact point enumerations are certified separately in rows 4--6a |

Two entries are deliberately not evidence and are labelled so rather than
promoted.  `sep_jacobian_audit.py` is a FALSIFIER: finite differences against
an interval enclosure can refute a Jacobian but never certify one, so it is a
regression guard, not a link in the chain.  `jacobian_fix_delta.json` records
which published numbers MOVED under the Jacobian fix; it justifies retracting
numbers, which is the safe direction, and the replacement numbers it names
are themselves in row 10 and await certificates.

Scope of every iv row (skill §1c): an interval certificate is about the
CONCRETE census map at the teachers it names, never the basis of a general
theorem, and it says nothing about any other teacher.

**The two compact domains are not the same domain, and the chain says so.**
Row 6b's `censusStripJ` is the enumeration's search rectangle written in the
labelled coordinates `(θ₀, θ₁)`: one period of `θ₁`, and a gap `θ₀ − θ₁` in
`[delta, dmax]` — a fundamental domain for the UNORDERED pairs, on which the
certificates of rows 4–6a are statements.  Row 8's `censusSquareJ` is the
LABELLED square, both angles in one period with a minimum gap, on which the
student swap acts and the factor two lives.  Transporting a count from the
strip to the square rests on the `2π`-periodicity of the angle map in each
student angle, which IS landed
(`CountConstancyJ/AnglePeriodJ.lean`: `censusAngleMapJ_add_two_pi_fst`,
`censusAngleMapJ_add_two_pi_snd`) — but the fundamental-domain bijection on
top of it is not, and no row claims it.

### The replayed witnesses

Each row is one iv certificate, replayed from its packed bitstring and the
spec alone: the leaves are re-derived, the boxes regenerated, every verdict
recomputed, and the solution enclosures re-derived for the boundary margin and
for the dedup.  Every row below reports zero failed leaves, zero undecided
leaves, and certified pairwise disjointness of the classes.  `families` is the
count of unordered separated families — the classes left after the fold and
`2π` merges, which are accepted only on a Krawczyk `unique` verdict for the
inflated joint hull — and does NOT include the exact fit, which is identified
by its own lattice test and counted separately.

| face | `(beta, y)` | families | leaves | boundary margin |
|---|---|---|---|---|
| F1 | `(0.75, 0.4)` | 4 | 254 893 | 1.37e-01 |
| F2 | `(0.3, −0.02)` | 3 | 416 339 | 1.37e-01 |
| F2 | `(0.75, −0.02)` | 3 | 171 381 | 1.37e-01 |
| F2 | `(1.5, −0.04)` | 3 | 73 458 | 1.37e-01 |
| F2 | `(2.2, −0.06)` | 3 | 56 443 | 1.37e-01 |
| F4 | `(0.3, −0.5)` | 1 | 84 199 | 1.30e-02 |
| F4 | `(0.75, −0.4)` | 1 | 59 233 | 5.83e-03 |
| F4 | `(1.5, −0.6)` | 1 | 46 281 | 1.37e-01 |
| F4 | `(2.0, −0.3)` | 1 | 31 817 | 1.37e-01 |
| F5 | `(2.4, −0.5)` | 3 | 26 301 | 5.00e-02 |
| F5 | `(2.6, −0.5)` | 3 | 26 524 | 5.00e-02 |
| F5 | `(2.8, −0.5)` | 3 | 28 783 | 5.00e-02 |
| F5 | `(3.0, −0.42)` | 3 | 60 088 | 1.37e-01 |
| F5 | `(3.0, −0.4)` | 1 | 60 890 | 1.37e-01 |
| F5 | `(3.0, −0.5)` | 3 | 58 421 | 5.00e-02 |
| F5 | `(3.0, −0.58)` | 3 | 60 019 | 1.37e-01 |
| F5 | `(3.0, −0.60)` | 1 | 60 860 | 1.37e-01 |
| F3 | `(3.08, −0.08)` | 3 | 273 978 | 2.65e-02 |
| F6 | `(0.75, −0.98)` | 3 | 166 091 | 1.37e-01 |
| F6 | `(1.5, −0.96)` | 3 | 72 195 | 1.37e-01 |
| F6 | `(2.2, −0.94)` | 3 | 53 570 | 1.37e-01 |
| F7 | `(3.08, −0.92)` | 3 | 226 693 | 2.65e-02 |

### F5 point enumerations have different total family counts

The F5 rows are the point: the certificates are valid and exhaustive at their
teachers, and they give different complete counts — three families at
`y = −1/2` and one at `(3.0, −0.4)` — under the same proposed F5
label. Thus the label cannot be treated as one total-count region. The selected
census does not change at these points because the extra families are saddles.

The drawing proposes a subdivision into a middle band and two outer strips,
but the point certificates do not prove global wall curves, their complement
components, or count constancy on those components. The conditional theorem
`censusAngleMapJ_familyCount_eq_witness_censusSquare` would transport a
count on a preconnected region only after its boundary-free and fold-free
premises were discharged there.

**At `beta = 3.0`, two pointwise count changes are bracketed**, each by a pair of
complete point certificates with different counts:

| candidate transition | below | above | point-count bracket |
|---|---|---|---|
| upper | 3 families at `y = −0.42` | 1 family at `y = −0.40` | at least one count-changing event occurs for `y ∈ (−0.42, −0.40)` |
| lower | 1 family at `y = −0.60` | 3 families at `y = −0.58` | at least one count-changing event occurs for `y ∈ (−0.60, −0.58)` |

with three families at the three enumerated heights `y = −0.58`, `−0.5`,
and `−0.42`, and one at the two enumerated outside heights. Every one of
those certificates is
valid and complete with zero undecided leaves, so the count really does change
inside each bracket; neither the number nor the type of intervening events is
identified. This does not enclose a wall curve. The conditional Lean
side of a possible subdivision is
`censusAngleMapJ_familyCount_eq_witness_on_censusBox`, which takes a
`beta`-interval times a `y`-interval and propagates a witness count once its
boundary-free and fold-free hypotheses are supplied. The complete point
enumerations alone do not supply those regional hypotheses. The
`beta`-dependence of the two candidate transitions — the drawing's band widens with
`beta` — is prospecting (`f5_count_band.json`) until the same bracketing is
done at a second `beta`.

Every certificate is in the table: **twenty-two witnesses,
twenty-two replays, all VALID AND COMPLETE**. A row entered only after its
replay exited `0`; the table was regenerated from the replay logs by
`fill_witness_table.py`, never typed.

## What is NOT closed

* **the separated locator's completeness** -- CLOSED AT ALL TWENTY-TWO NAMED
  WITNESSES, open as a complete-global-label statement. The sign-law-pruned
  enumeration (`face_bb_signlaw.py`)
  finishes with zero undecided area at every teacher in the replay table:
  every point of the `(th1, D)` domain is either excluded by a certificate
  or inside a Krawczyk box, so the certified family list there is complete —
  exactly, not "at least" — and the count agrees with the locator at every
  named witness (`witness_closure.json`). What remains open is transport from
  a witness to a whole proposed label, except on the separately certified
  68-tile rectangle. That transport needs the regional boundary and
  fold-freeness inputs of the count-constancy schema;
* **F5's two pointwise count transitions** -- the family count changes at
  `beta = 3.0` across two straddling pairs of complete certificates.
  What remains open is whether these samples belong to global wall curves and
  how any such curves depend on `beta`; the apparent widening in
  `f5_count_band.json` is prospecting. At the replayed points and continuation
  grid, the extra separated families are saddles and therefore do not change
  the selected census string. Constancy of that string across a complete
  global F5 cell is not claimed;
* **F1 is exact at its named teacher, not on the proposed label.** The packed
  replay at `(0.75,0.40)` exhausts the search domain with 254,893 leaves
  and four separated families. Earlier 1%-residual branch-and-bound
  measurements are historical prospecting, not the current status. Transport
  across the proposed F1 label remains open outside the separately certified
  atlas;
* **F3 and F7 are exact at their named teachers.** Both had been abandoned as
  permanently plateaued
  under the old Cramer Jacobian — F3 at `1.4e-2`, F7 at `1.2e-2`, unmoving
  after 14 000 s. At their named `beta = 3.08` witnesses the smallest gap is
  `0.0315`, so `delta = 0.02` would already cover every separated family.
  The shipped runs used the finer `delta = 0.005`; that finer collar is useful
  for continuing the column toward `beta = pi`, where the smallest gap drops
  below `0.02`, but it was not required for exactness at the named witnesses
  themselves (see the validity-domain note beside the spec). Retried at that
  finer floor with the division-free centered form, both finished
  with **zero undecided boxes**: F3 at `(3.08, −0.08)` in 4 956 s, 273 978
  leaves, 3 separated families; F7 at `(3.08, −0.92)` in 4 081 s, 226 693
  leaves, 3 separated families. Together with F1, F2, F4, F5, and F6, every
  proposed label now has at least one exhaustively enumerated teacher. Neither
  point needed a new idea — they needed the division-free derivative, the
  same repair that broke the mixed-sector wall;
* **row 7a's mass-free route: three holes and one cost wall**, all named and
  none of them an enclosure defect.  (i) The exact fit, where
  `generalJKernelTeacherSchurDetT` vanishes identically — now a Lean theorem,
  and it makes the un-restricted exclusion FALSE rather than open, so row 7b's
  mass-carrying route is not an alternative there but the only route.  (ii) The
  antipodal locus `θ₀ − θ₁ = π`, which holds five of the 55 certified
  non-exact-fit families EXACTLY, all on `y = −1/2`; they are outside
  `GeneralJSeparatedPairCritical` by that predicate's own non-opposition field,
  so the mass-free route cannot reach them at any precision.  (iii) The torque
  conditions, which the tree can only supply region-wise under both masses
  positive. The only surviving row-7a interval certificate is the map-only
  `row7a_smoke`; no named proposed-label rectangle has a valid Schur
  exclusion. Prospecting estimates about `1.5e9` boxes for one proposed
  face and gap band, and no band reaching below a student gap of about
  `0.4` appears coverable — there the determinant test needs box sides
  from `1e−5` down past `1e−8`;
* **the mixed collar samples** — located families certified individually;
  only `beta = 0.15` has a two-sided locate-then-certify transition
  pair, while the other eight samples provide one-sided first-success
  heights; no global wall is proved;
* **uncertified strips** hugging the strata (`y` below ~0.004–0.03,
  `beta`-dependent) and 19 declined antipodal-band points, all with
  `|beta − pi| <= 0.0240`. Every sampled point with
  `|beta − pi| >= 0.0376` certifies and agrees with the shipped
  classifier; certification inside the pocket is mixed;
* the columns `beta = 0`, `beta = pi` are exact degenerate cases, treated
  in closed form, not gaps.

The coincident stratum has NO unresolved set: the analytic closure removed
the six collars and the two never-terminated outer coverings of the first
version of this certificate.

**Why the gap floor `delta > 0` is not a convenience.**  The whole
coincidence diagonal `θ₀ = θ₁` lies inside the zero set of the angle map, at
every teacher: there the kernel takes its diagonal value `phiCosJ 0 = π`, so
the radial determinant vanishes, and the coupling `couplingHJ 0` vanishes
too, killing both equations
(`SeparatedCount/SwapPairingJ.lean`: `separatedAngleMapJ_eq_zero_of_diag`).
A counting domain that meets the diagonal therefore has an infinite zero set
and no count at all — so the enumeration's `D ≥ delta > 0` is what makes the
count-constancy schema applicable, not merely what makes the search
terminate.

## Validation

25 000 samples against the literal function sources of `website/widgets.js`
at commit `2201dd0`: kernel atoms ≤ 8.9e-16, Wronskians ≤ 2.7e-15, `Dsep`
1.8e-14, `y` 4.4e-16 — all 1–3 ulp, and **bit-identical** when node's own
`sin`/`cos` are fed into the Python formulas (the residue is libm's,
measured at 1.1e-16).  `h' = slopeAtomJ` by 300-digit differences: 5.4e-51.
The analytic `Dsep` gradient against differences: 3.3e-39.

## Reproducing

Requires `mpmath`, `sympy`, `node`.  From this directory:

```
node ref_widgetsJ.js               # reference samples from the literal widgets.js
python3 validate.py                # widgets.js agreement, reduction, identities
python3 certify_analytic.py        # the coincident closure  -> analytic_counts.json
python3 spot_check.py              # outer-interval corroboration
python3 conditioning_measure.py    # the sweep's conditioning (slow)
python3 wall_pos.py 0 1            # the mixed collar scan (slow)
python3 faces.py 0 1               # the 28 face witnesses (slow)
python3 sep_jacobian_audit.py      # REGRESSION GUARD for the Krawczyk Jacobian
python3 jacobian_fix_delta.py all  # item-by-item delta of the Jacobian fix
python3 faces_rebuild_fixed.py     # post-fix witness table -> faces_fixed.json
python3 band_rescan_fixed.py 0 1   # post-fix band          -> bandfix2_*.json
python3 f5_count_band.py           # prospect F5's internal count band and candidate transitions
python3 face_bb_signlaw.py <beta> <y> <budget> 1 <tag> - - <delta> hybrid
                                   # the sign-law-pruned enumeration (slow)
python3 witness_closure_table.py   # per-witness closure -> witness_closure.json
python3 cont_gate.py               # FALSIFIER for the angle-map Jacobian determinant
python3 cont_cert.py F4smoke       # the continuation round trip, a few seconds
python3 cont_check.py cont_F4smoke.json
python3 cont_cert.py F4            # the F4 continuation certificate  (~2 min)
python3 cont_check.py cont_F4.json # replay, exit 0 = valid and complete
for f in F1a F1b F2 F3 F5a F5b F6 F7; do   # the other eight rectangles
  python3 cont_cert.py  $f
  python3 cont_check.py cont_$f.json
done                               # ~85 min of emit in total, one core
python3 kernel_schur_mv_gate.py    # GATE for the mean-value determinant form
python3 validate_rowform.py        # GATE for the division-free row form
python3 row7a_band.py              # required box side per family, beta thin AND swept
python3 row7a_stagewidth.py        # where the small-gap enclosure is lost, by stage
python3 row7a_stagewidth.py 160 all  # the dependency factor at every family
python3 row7a_rowgate.py           # GATE: the MAP verdict tests the census angle map
python3 row7a_check.py row7a_smoke.json       # replay, exit 0 = valid and complete
python3 row7a_check.py row7a_safesmoke.json   # expected regression failure
```

`cont_gate.py` must be re-run after ANY edit to `census_cert.sep_res_rowform_jac`:
it is the only check that can see an invalid Jacobian-determinant enclosure,
and that enclosure is what every continuation leaf's fold-freeness claim rests
on.

`sep_jacobian_audit.py` must be re-run after ANY edit to `census_cert.sep_res`:
it is the only check that can see an invalid derivative enclosure, which is
what CERTIFICATE 9.5 records.

`certificate.json` is the committed record of the original box-covering run;
its per-run input shards were not retained, so `assemble.py`/`finalize.py`/
`fill_md.py` do not re-run from the committed tree (Files section of
`CERTIFICATE.md`).  `ref_samplesJ.txt` is regenerated by `ref_widgetsJ.js`
rather than committed.
