# Selected-minimum DFS certificate

`minimum_dfs.py` is the numerical consumer of the Lean selected-zero theorem
in `CountConstancyJ/MinimumSelectorJ`.  It certifies only the roots that can be
strict local minima:

```text
censusAngleMapJ = 0,  CramerT00 >= 0,  CramerDetT >= 0.
```

Every replayed leaf either excludes an angle-map zero, proves that the first
pivot or determinant is negative (so the root is unselected), or proves the
determinant strictly positive.  A saddle fold therefore does not become a
minimum-map wall.

## Exact formal bridges

The certificate uses the following named Lean facts.

- `censusMinimumSelectedZeroCountJ_eq_on_face` says the selected-zero count is
  constant when selected roots avoid coincidence/opposition and their Cramer
  determinant is nonzero.
- `censusMinimumSelectedZeroCountJ_eq_witness_strip_of_certificates` is the
  exact end-to-end form: a positive gap floor and `dmax < pi` derive the two
  lattice exclusions, while the regularity tree, four boundary trees, and
  witness count supply its remaining hypotheses.
- `censusMinimumSelectedZeroCountJ_eq_witness_strip_of_relativeCertificates`
  is the production-chart specialization.  `RelativeChart.lean` normalizes
  `u=theta1-beta` by an integer full turn and proves the angle map, pivot, and
  determinant invariant, so a relative full-period tree discharges regularity
  on the fixed physical strip.
- `MAP_AUX` additionally uses the proved positive-mass straddle and
  mixed-mass same-side laws after both Cramer mass numerators exclude zero;
  the same code also records parameter-uniform Krawczyk exclusions.
- `generalJCramerSchurT00_factor` and
  `generalJCramerSchurDetT_factor` identify the factored interval expressions
  with the formal selectors exactly.
- `censusNormalizedAngleMapJ_eq_zero_iff` identifies the small-gap
  `D^-3` Taylor form with the original angle-map zero set for `D != 0`.
- `censusMinimumSelectedZeroJ_isStrictLocalMin` realizes a selected regular
  zero as a strict local minimum of the original noncentered loss.

All named statements are axiom-audited.  Their allowed dependency list is
`[propext, Classical.choice, Quot.sound]`.

## Certificate format and replay

The domain coordinates are `(beta,y,s,D)` with
`theta1=s, theta0=s+D`, or the teacher-aligned chart `(beta,y,u,D)` with
`theta1=beta+u, theta0=beta+u+D`.  The latter removes dependency along the
moving kernel kink `theta1-beta=0`.

The artifact stores only a parameter specification and a preorder bitstream.
The checker regenerates every root and child, recomputes every interval
verdict by its named proof method, checks total covered volume, and rejects any
`UNDECIDED` leaf.  Replaying the named method rather than the search's current
priority order keeps an older tree valid when a sharper exclusion is inserted
ahead of it.  The open separated domain is restricted to `0 < D < pi`;
coincidence and opposition are separate strata.

Boundary containment is a separate replay obligation.  The regularity sweep
may legitimately accept a `POS_DET` leaf, but a physical search face must
contain no selected root at all.  `minimum_boundary_dfs.py` therefore works in
absolute coordinates and accepts only map exclusions or negative selectors on
the four faces `seam_lo`, `seam_hi`, `gap_lo`, and `gap_hi`.
`MinimumSelectorJ/Boundary.lean` turns those four exclusions into the strict
inequalities consumed by the strip count theorem.

Run the falsification gate and replay an artifact with:

```bash
python3 minimum_dfs_gate.py
python3 minimum_dfs_check.py minimum_dfs_f5fold_smoke_p000.json
python3 minimum_boundary_check.py minimum_boundary_<tag>_<face>_p000.json
python3 minimum_witness_check.py minimum_witness_f4probe_census.json
python3 minimum_atlas_check.py minimum_atlas_f4.json --structure-only
```

For the full compact F4 rectangle, `minimum_atlas_sweep.py` schedules the
`0.1`-wide connected grid (with two `0.05` edge rows), resumes from discovered
artifacts, and records a tile only after replaying its regularity forest and
all four boundary forests.  The journal is not a certificate; after every tile
passes, assemble and independently replay the ordinary atlas manifest.  The
atlas checker also subdivides the declared target at every box endpoint and
rejects any uncovered parameter cell:

```bash
python3 minimum_atlas_sweep.py --workers 20
python3 minimum_atlas_sweep.py --write-atlas minimum_atlas_f4_full.json
python3 minimum_atlas_check.py minimum_atlas_f4_full.json --workers 20 \
    --dfs-checker minimum_dfs_check_cached.py
```

The cached replay entry point starts a fresh process, installs the same
bounded exact atom cache used by `minimum_dfs_cached.py`, and then calls the
unchanged `minimum_dfs_check.main`.  The cache key includes the interval
precision and mpmath's exact directed endpoints; it changes no bit reader,
box, split, target proof method, coverage test, or acceptance condition.
`minimum_dfs_check_cached_test.py` runs both fresh entry points on the same
valid forest and on the same invalid forest and requires identical output and
exit status.  Use `minimum_dfs_check.py` directly when an uncached replay is
desired.

One exceptionally large ordinary forest can instead be replayed exactly over
a deterministic tree frontier:

```bash
python3 minimum_dfs_check_parallel.py minimum_dfs_<tag>_p061.json \
    --depth 8 --workers 8
```

The coordinator strictly parses the entire preorder forest and stored
statistics, stopping at depth 8 or at an earlier leaf.  Those derived paths
are a prefix-free cover; no path list is accepted from the artifact.  Every
spawned worker reopens the SHA-256-identified file, reconstructs the same
paths and boxes with `minimum_dfs.child`, and calls the same targeted
`minimum_dfs.verdict` proof-method check as the serial checker.  The
coordinator rejects missing or duplicate assignments and recomposes all bits,
nodes, verdict counts, leaf volumes, and root volume before accepting.  The
artifact remains an ordinary `minimum_dfs` forest and the serial checker can
always replay it.  Pass `--no-cache` to use this decomposition without exact
atom memoization.  `minimum_dfs_check_parallel_test.py` compares fresh serial,
cached-parallel, and uncached-parallel replays and adversarially checks
truncation, extra bits, padding, false claims, undecided leaves, stored stats,
prefix assignments, and artifact hashes.

Atlas `--workers` already parallelizes distinct witness, regularity, and
boundary artifacts.  When a few large regularity forests are the long pole,
the parallel checker can be nested deliberately; keep the product of outer
and inner workers within the host CPU budget:

```bash
MINIMUM_REPLAY_WORKERS=4 python3 minimum_atlas_check.py \
    minimum_atlas_f4_full.json --workers 4 \
    --dfs-checker minimum_dfs_check_parallel.py
```

## F4 rectangle closed with count one

The artifact `minimum_dfs_f4probe_p000.json` covers

```text
beta in [0.7,0.8], y in [-0.45,-0.35],
u over one full period, D in [0.001,3.13].
```

Its independent replay consumes all 473,671 bits and all 94,747 leaves:
83,861 angle-map exclusions, 100 negative-selector leaves, 10,786
positive-determinant leaves, and zero undecided leaves.  Hence every selected
root in this rectangle has nonzero determinant.  The relative-chart sweep is
transported to the physical chart by the angle-map period theorem and the
matching selector-period theorems in `MinimumSelectorJ/Period.lean`.

The four absolute-coordinate artifacts
`minimum_boundary_f4probe_{seam_lo,seam_hi,gap_lo,gap_hi}_p000.json` replay
1,403 leaves and prove that the selected set misses every physical face of
`censusStripJ 0.137 0.001 3.13`.  The bridge in
`MinimumSelectorJ/Boundary.lean` converts those exclusions into the strict
interior hypothesis of the selected-count theorem.

Finally, `minimum_witness_f4probe.json` certifies the count at
`(beta,y)=(0.75,-0.4)`.  A 34-leaf normalized collar excludes all roots for
`D in [0.001,0.02]`; the existing 59,233-leaf complete enumeration supplies
the rest.  After removing its one fold copy above `D=3.13`, two certified
root classes remain.  Direct interval evaluation of the formal selectors
gives

```text
D = 0.750000000:  T00 > 0, det > 0  (selected exact fit)
D = 3.097422789:  T00 > 0, det < 0  (unselected saddle).
```

Therefore `censusMinimumSelectedZeroSetJ` has cardinality one throughout this
named F4 rectangle.  This is a closed rectangle of the minimum map, not yet a
certificate for the whole nonrectangular F4 face or for its adjacent strata.

Reproduce the three replay layers from this directory with:

```bash
python3 minimum_dfs_check.py minimum_dfs_f4probe_p000.json
for f in minimum_boundary_f4probe_*.json; do
  python3 minimum_boundary_check.py "$f" || exit 1
done
python3 minimum_witness_check.py minimum_witness_f4probe.json
```

## Selected-regularity smoke near the F5 equal-mass line

The eight `minimum_dfs_f5_selected_relative_pNNN.json` shards cover

```text
beta in [2.59,2.61], y in [-0.51,-0.49],
u over one full period, D in [0.02,3.13].
```

All shards replay with zero undecided leaves. They certify selected regularity
on exactly those product boxes: a box is discharged by an angle-map exclusion,
a negative selector, or a nonzero determinant. The smaller
`minimum_dfs_f5fold_smoke_p000.json` is another selected-regularity smoke inside
that same near-equal-mass range.

These artifacts do not certify a fold, a selected count, or a neighbourhood
closure. In particular, the accepted continuation point at `beta = 2.60`,
`y = -0.4580` lies outside their teacher range. A regional count would also
need the four physical-boundary certificates and a complete witness count used
by the end-to-end bridge above; those layers are absent here.

## Census scope: counting EVERY zero, not only the minima

`spec["scope"]` selects which zeros the certificate controls.

* `"selected"` (the default, and what every artifact without the key means):
  only zeros that can be local minima.  `NEG_T00` discharges a box by proving
  its first Cramer-cleared pivot strictly negative, so any zero inside is an
  unselected saddle.  A saddle-saddle fold is therefore not a wall.
* `"census"`: every zero of the angle map.  Being a saddle is no longer an
  excuse, so `NEG_T00` is not an admissible discharge — only an angle-map
  exclusion or a definite determinant sign closes a box.  The boundary DFS
  narrows the same way, to map exclusions only, because a negative pivot says a
  zero is a saddle and not that the face is free of it.

Census scope exists because it is what the face-count schema consumes.
`censusAngleMapJ_zeroCount_eq_witness_strip` needs fold-freeness at every zero,
which `separatedAngleJacDetJ_ne_zero_of_schurDetT_ne_zero` reduces to
`det T != 0`; `FaceClosureJ/CramerFoldFreeJ.lean` supplies both that and the
all-zero boundary input from the replayed artifacts, as
`censusAngleMapJ_zeroCount_eq_witness_strip_of_certificates`.  A constant count
across a region, anchored by one COMPLETE enumeration at a witness, is what
discharges the separated locator's completeness obligation there: a family the
locator never seeded would make the count exceed the certified value.

A census forest is also a valid selected certificate; a selected forest is NOT
a valid census one.  `tile_id` therefore puts the scope in the artifact tag and
`discover(scope)` filters on it, so the weaker artifact can never stand in for
the stronger.

### The F4 probe rectangle, closed at census strictness

`beta in [0.7,0.8]`, `y in [-0.45,-0.35]`, strip `seam=0.137 D=[0.001,3.13]`.

```text
minimum_dfs_f4probe_census_merged.json     64 shards, 473,726 bits, 94,758 leaves
   MAP_CENTERED 70,746   POS_DET 10,786   MAP_COLLAR 7,094   MAP_AUX 5,414
   MAP_PLAIN       627   NEG_DET     91   NEG_T00        0   UNDECIDED   0
minimum_boundary_f4probe_census_{seam_lo,seam_hi,gap_lo,gap_hi}_p000.json
   1,541 leaves, map exclusions only
```

`NEG_T00 = 0` is the content: no zero was set aside for being a saddle.  All
five artifacts replay `exit=0`, leaf volume matching root volume to 16 digits.

Reproduce with:

```bash
python3 minimum_dfs_check.py minimum_dfs_f4probe_census_merged.json
for f in minimum_boundary_f4probe_census_*_p000.json; do
  python3 minimum_boundary_check.py "$f" || exit 1
done
python3 minimum_witness_check.py minimum_witness_f4probe_census.json
```

The witness replay reports 2 zeros in the strip at `(0.75,-0.4)`: the exact fit
at `D = 0.750000000` (`T00 > 0`, `det > 0`, selected) and an unselected saddle
at `D = 3.097422789` (`T00 > 0`, `det < 0`), both with `det != 0`.

**Cost.** Against the selected forest on the identical rectangle (specs equal
but for the scope key): 94,758 leaves versus 94,747, i.e. **0.012 % more**.  The
nine `NEG_T00` boxes subdivided into twenty angle-map exclusions and nothing
else moved.  Across seven earlier forests `NEG_T00` was reached on 0.048 % of
leaves and 53.6 % of those already satisfied a determinant sign on the same box.

**Scope.** This fixes the ZERO count on one rectangle.  Reading it as a FAMILY
count goes through the census lane's `#zeros = 2*(N+1)` identification, which
the schema explicitly does not assert.  On this strip the measured relation is
`#zeros = N_separated + 1`, because `dmax = 3.13 < pi` excludes the `2pi-D`
copy — the witness enumeration shows exactly that, 3 raw enclosures with 1 fold
copy outside.  Measured at one teacher; not a landed identification.

Grow a CONNECTED census region (the closure quantifies over a preconnected
region, so scattered census tiles buy nothing):

```bash
python3 minimum_atlas_sweep.py --workers 8 --sharded --shard-size 1 \
    --shard-budget 300000 --scope census --beta-min 0.5 --y-min -0.45 \
    --dfs-script minimum_dfs_cached.py \
    --dfs-checker minimum_dfs_check_cached.py
python3 minimum_atlas_sweep.py --scope census --beta-min 0.5 --y-min -0.45 \
    --write-atlas minimum_atlas_f4_census_safe.json
python3 minimum_atlas_check.py minimum_atlas_f4_census_safe.json --workers 8 \
    --dfs-checker minimum_dfs_check_cached.py
```

The sweep generator and regularity checker are separate allowlisted choices.
The cached checker installs the exact atom cache and then calls the unchanged
ordinary checker, so it changes neither the bit reader nor any acceptance
condition; boundary generation and replay keep their ordinary entry points.
For an additional independent audit, replay the assembled manifest again with
`--dfs-checker minimum_dfs_check.py`.

Two spare cores can precompute far-queued root cells without creating a second
journal writer:

```bash
python3 minimum_atlas_aux_root.py --scheduler-pid PID \
    --expected-workers 8 --scope census \
    --root 2.1 2.2 -0.25 -0.15 48 \
    --root 2.0 2.1 -0.15 -0.1 48
```

This is deliberately not another sweep.  Each root is first written under a
unique staging tag.  The driver then stops only the scheduler parent, requires
all eight direct children to be recognized, refuses any active target tile,
rechecks the exact spec/range/scope and zero `UNDECIDED` count, and publishes
with an atomic no-clobber hard link.  The scheduler is resumed in a `finally`
block even when a guard refuses publication.  Auxiliary roots never touch the
journal, and active tiles are forbidden because their missing ranges were
snapshotted before the auxiliary file existed.

For sustained spare-capacity use, `minimum_atlas_aux_loop.py` repeatedly sends
exactly two distinct far-queue tiles through that same guarded publisher.  It
chooses each tile's current first missing root and stops on host-load pressure,
when the live sweep frontier gets within the configured queue gap, or when any
guard refuses a batch.  It also never reads or writes the journal.
For a coordinated scheduler cutover, pass `--stop-file PATH` and create that
sentinel: the loop samples it only after a batch has been fully published and
the scheduler resumed, then exits before staging another root pair.

When the two target rows have strongly asymmetric root costs,
`minimum_atlas_aux_lanes.py` avoids leaving the faster slot idle.  It publishes
completed roots one at a time through the same serialized guard and refills
that lane only after independently rechecking `load1 < --max-load`, exact main
worker accounting, the frontier gap, active targets, and no-clobber.  Its stop
sentinel drains and publishes existing lanes without starting replacements.

Both cutoffs are part of the claim: together they select the closed 68-tile
rectangle `beta in [0.5,2.2]`, `y in [-0.45,-0.1]`.  The larger 84-tile
rectangle starting at `beta=0.1` is not admissible: the measured `D=dmax`
antipodal exit curve crosses its bottom row for `0.3 <= beta < 0.5`, while the
curve below `beta=0.3` is unmeasured rather than absent.  The scope-aware atlas
checker requires census artifacts on every box and anchors their common full
zero count at the two-zero witness; it cannot consume a selected-only forest.

## Closure of the 68-tile census rectangle (measured 2026-08-17)

The sweep over `beta in [0.5,2.2]`, `y in [-0.45,-0.1]` is COMPLETE.  All 68
tiles are journaled `passed` in `minimum_atlas_f4_sweep_census_state.json`,
which records a tile only after a FRESH PROCESS has replayed its merged
regularity forest and all four physical boundary forests.  The manifest
`minimum_atlas_f4_census_safe.json` assembles those 68 boxes and passes

```bash
python3 minimum_atlas_check.py minimum_atlas_f4_census_safe.json --structure-only
```

with `ATLAS STRUCTURE VALID (REPLAY SKIPPED): manifest claims angle-map zero
count = 2; the count was not independently replayed`, a single connected chain
from the base box `f4probe_census` through all 68 boxes, and
`target: beta=[0.5,2.2], y=[-0.45,-0.1], fully covered`. The structure-only
call checks the manifest graph and geometric coverage; the count evidence is
the separately recorded fresh per-tile replays described above, not the
structure-only verdict itself.

The concrete tile constants and their replay verdicts remain external data:
the Lean atlas theorem is conditional and generic, and does not contain these
68 tiles. The identification of the Python checker's interval expressions
with the Lean constants `censusAngleMapJ`, `censusCramerSchurT00J`, and
`censusCramerSchurDetTJ` is human-checked rather than a Lean theorem. Thus the
regional conclusion composes the external replay and this stated semantic
correspondence with the conditional Lean propagation theorem.

Measured totals over the 68 merged census forests:

```
preorder nodes      206,853,513
leaves               41,371,573
  MAP_CENTERED       36,356,266     angle-map exclusion, centered branch
  POS_DET             3,363,051     definite positive determinant
  MAP_COLLAR          1,078,635     angle-map exclusion, collar branch
  MAP_AUX               478,915     angle-map exclusion, auxiliary branch
  MAP_PLAIN              57,081     angle-map exclusion, plain branch
  NEG_DET                37,625     definite negative determinant
  NEG_T00                     0     inadmissible in census scope, and absent
  UNDECIDED                   0
unresolved volume             0.0
replay wall time      460,454 s (127.9 h), slowest tile 45,875 s
artifact volume          34.5 MB across 68 merged forests
```

`NEG_T00 = 0` is a consequence of the scope, not a coincidence: a census
forest may close a box only by an angle-map exclusion or a definite
determinant sign, so every one of the 41,371,573 leaves is one of those two
kinds, and `UNDECIDED = 0` with zero unresolved volume means the closed
rectangle is exhausted.

**What was checked, exactly (nothing more).**  Every artifact under the
manifest was replayed ONCE, by a fresh-process checker: per tile the ordinary,
cached, parallel, or dynamic-parallel entry point, as recorded in the journal.
The assembled manifest was then checked structurally.  The optional SECOND,
redundant replay of all 68 forests from the manifest,

```bash
python3 minimum_atlas_check.py minimum_atlas_f4_census_safe.json --workers 8 \
    --dfs-checker minimum_dfs_check.py
```

was NOT run for this instance; it costs another ~128 core-hours and re-derives
what the per-tile replays already established.  Two tiles present in the
journal lie OUTSIDE the declared rectangle (`y in [-0.55,-0.45]`,
`beta in [0.7,0.9]`) and are excluded from the manifest and from the claim.
