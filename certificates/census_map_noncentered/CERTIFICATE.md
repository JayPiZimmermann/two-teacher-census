# Evidence ledger for the NONCENTERED (plain-ReLU) two-student census map

**Status date: 2026-08-24.** This file is the authoritative claim ledger for
the noncentered census certificate. `README.md` is the reproduction and
engineering guide; generated JSON objects point back here for their scope.

**Object.** Analytic coincident-boundary results, sampled candidate-transition windows,
pointwise and regional certificates, and the proposed incidence graph used to
draw the noncentered census map, in the coordinates
`(x = beta in [0, 2pi), y in (-1, 1])` of `website/widgets.js`.  This document
does not certify that the proposed graph is the complete global arrangement;
section 7 states that boundary explicitly.  Companion to
`../census_map_centered/`; the harness, validation method, and output format
are shared.

**Definitions.** The teacher gap is `beta in [0, 2pi)` and the mass pair is
`(s0, s1) = (sin psi, cos psi)`, `psi = (y+1) pi/2` (`massesAt`); negating
both masses is an exact isometry, so `y = -1` and `y = +1` are the same line
and `y` lives on a circle of circumference 2 (`ratioCoord` folds by it).
The kernel is the Lean `phiCosJ y = reluPhi (cos y)` with
`reluPhi rho = orthantPhi rho + (pi/2) rho`, and its coupling
`couplingHJ y = couplingH y + (pi/2) sin y`
(`Planar/SignedAngle/PlanarKernelJ.lean`); unlike the centered kernel both
are only `2pi`-periodic.  The closed branch forms used by `widgets.js`
(lines quoted are functions, matched by identifier, not line number) and
mirrored in `noncentered.py` are

```
  phiJ(t) = (pi - x) cos x + sin x ,  x = mod(t, 2pi) folded by x > pi -> 2pi - x
  HJ(t)   = |pi - x| sin x         ,  x = mod(t, 2pi)          [= couplingHJ]
  dHJ(t)  = HJ'(t) = phiJ(t) - 2|sin t| = slopeAtomJ
  kappa = pi ,  period = 2pi
```

and agree with the Lean definitions for every `t` (on `x in [0, pi]`,
`arcsin(cos x) = pi/2 - x` turns `reluPhi(cos x)` into `(pi - x) cos x +
sin x`; evenness and `2pi`-periodicity carry the rest).  The three mass-free
determinants are built from these atoms exactly as in the centered
certificate:

```
  Wtau(beta,t) = h(t) sA(t-beta) - h(t-beta) sA(t)          torque Wronskian
  Wpot(beta,t) = phi(t) h(t-beta) - phi(t-beta) h(t)        potential Wronskian
  Wwgt(beta,t) = |sin t| h(t-beta) - |sin(t-beta)| h(t)     weight Wronskian
```

with `h = HJ`, `phi = phiJ`, `sA = dHJ`.  A degenerate torque root `t` whose
kernel vector `(-h(t-beta), h(t))` is nonzero maps to the `y` coordinate
through that mass direction, folded exactly as `ratioCoord`.  At a common-zero
pair the vector and ratio coordinate are undefined; degenerate columns are
handled separately.

**Precision.** All reported enclosures come from `mpmath.iv` (interval
arithmetic with directed rounding), mpmath 1.3.0, at **200 bits** for the
coincident stratum, **160 bits** for the per-teacher census, **120 bits** for
the analytic-closure inequality, **90 bits** for the two-dimensional
mass-free sweep.  `float64` is used only in the `widgets.js` mirror of
section 1 and in *locators* that guess where a solution is; no reported
enclosure passes through it.

---

## 0. Headline

| evidence object | outcome | scope |
|---|---|---|
| proposed incidence graph on `beta in [0, pi]` | `V = 11`, `E = 18`, `F = 7`; conditional Euler bookkeeping `V-E+F=0` | assumes the interpolated candidate-transition arcs have the proposed incidences and that no additional wall exists |
| proposed full-torus doubling | `V = 14`, `E = 28`, `F = 14` | the same conditional graph, not a global arrangement certificate |
| representative pointwise censuses and collar row | five exact values across the seven replayed F1--F7 representatives; one additional collar pattern whose listed family is certified to exist | collar lists are not exhaustive; no set of global region values is proved |
| exhaustive separated enumerations | twenty-two named teachers, 2,368,457 leaves, zero failed or undecided leaves | exact pointwise lists only |
| regional atlas | 68 tiles, full zero count `2`, selected count `1` | `beta in [0.5,2.2]`, `y in [-0.45,-0.1]`, and `censusStripJ 0.137 0.001 3.13`; trust boundary in 9.7 |
| components drawn by commit `2201dd0` at depth 8 | 49, of which 19 are below the naming threshold | finite-resolution renderer/classifier measurement |

The shipped classifier of commit `2201dd0` -- whose separated locator is the
well-conditioned enumeration this certificate introduced (section 4.3),
adopted into `widgets.js` -- agrees with the returned family lists at all
**137 points where every returned family passes its interval check** (section 8; the pre-Jacobian-fix figure
146 of 146 is superseded by 9.5) and at all 28 witnesses in
`faces_fixed.json` (two relocated by 9.5).  The exhaustive replay makes the
lists exact at the twenty-two teachers in its ledger; agreement elsewhere is
with the returned, individually certified families and does not exclude an
unseeded one.  The predecessor classifier at commit `f430fab` agreed at 94 of
146 points and drew 1,868 sub-threshold fragments; section 8 records both
measured states against their commits.

---

## 1. Validation against `website/widgets.js`

Reference values are produced by `ref_widgetsJ.js`: `node` evaluating the
**literal function sources** of `website/widgets.js` (`PI`, `mod`, `phiJ`,
`HJ`, `dHJ`, `slopeAtom`, `Pot`, `Trq`, `ratioCoord`), extracted by
identifier with a guard asserting each extracted function still contains the
expected expression, and `eval`-ed verbatim.  Last re-run against
`widgets.js` at commit `2201dd0`.  Sample count: **25000** quasi-random
`(t, beta, t1)` triples, `t in [-2pi, 4pi]`, `beta, t1 in [0, 2pi]`.

| quantity | max abs deviation, Python mirror vs widgets.js |
|---|---|
| `phiJ` | 8.882e-16 |
| `HJ` | 4.441e-16 |
| `dHJ` | 8.882e-16 |
| `slopeAtom` | 8.882e-16 |
| `Wtau` | 2.665e-15 |
| `Wpot` | 2.665e-15 |
| `Wwgt` | 1.332e-15 |
| `Dsep` (separated determinant, division form) | 1.776e-14 |
| `y = ratioCoord(-h(t-beta), h(t))` | 4.441e-16 |

These deviations are **1-3 ulp and are entirely attributable to libm**
(`ref_libm.js` + section (a2) of `validate.py`): V8's `Math.sin`/`Math.cos`
and glibc's differ by up to `1.110e-16` (6000 arguments), and feeding node's
*own* `sin`/`cos` into the Python formulas reproduces every kernel quantity
**bit-identically (max deviation 0.0)** — the reimplementation is
formula-identical to the widget.

`h'(t) = slopeAtomJ(t)` by central differences (mpmath prec 300, step
`1e-25`, kink lattice `t in pi Z` avoided by `1e-3`, 799 samples): max
absolute error **5.4027e-51**; it is an exact identity of the branch forms
(section 2.2).  The analytic gradient of the separated determinant against
central differences (`mpmath.iv` prec 260, step `1e-20`, 175 samples): max
error **3.34e-39**.

Since commit `2201dd0` the shipped separated locator (`sepMassesJ`,
`sepBalanceJ`, `separatedScanJ`) uses the same radial-row elimination as this
certificate's section 4.3 — `massDet = pi^2 - phiJ(D)^2`, weights
`c0 = (kap P(th0) - phi(D) P(th1))/massDet` (and symmetrically), residual
`(c1 h(D) + A(th0), c0 h(D) - A(th1))` — formula for formula against
`census_cert.sep_res`.

---

## 2. An exact reduction that removes the kinks

### 2.1 The kink lattice is FOUR points, not two

`phiJ` and `HJ` are both `C^1` everywhere (`phiJ' = -HJ`,
`HJ' = slopeAtomJ`, both continuous), and `HJ` is `2pi`-periodic while
`|sin|` stays `pi`-periodic.  The atoms that carry a kink are
`slopeAtomJ = phiJ - 2|sin|` and `|sin|` itself, and they kink at
`t in pi Z`.  The determinants use them at `t` and at `t - beta`, so the
kink lattice in `t` over one `2pi` period is

```
      t = 0 ,   t = beta ,   t = pi ,   t = beta + pi        (mod 2pi)
```

**four** points, not the centered lane's two.  The extra pair — the kink at
`t = beta + pi`, with no centered counterpart — invalidates any interval
derivative test across it.  Every certified derivative test below is run on
a piece whose closure meets no kink in its interior, with the branch indices
supplied explicitly rather than inferred.

### 2.2 The closed branch forms, and the reduction

Index the half-branches by `n` with `t in [n pi, (n+1) pi]`, and put

```
   p = pi + 2pi*floor(n/2) - t          (so x = mod(t, 2pi) = pi - p)
   e = +1 if n even else -1             (e = sign(p), p in [0,pi] / [-pi,0])
```

Then the closed branch forms of `widgets.js` are exactly

```
   h(t)      = |p| sin p
   phi(t)    = -|p| cos p + e sin p  = e (sin p - p cos p)
   sA(t)     = -|p| cos p - e sin p  = phi(t) - 2|sin t|
   |sin t|   = e sin p
```

which makes `sA = phi - 2|sin|` an identity of the branch forms rather than
a numerical coincidence, and gives `h' = sA`, `phi' = -h`,
`sA' = 2 e cos p - |p| sin p`, `(|sin t|)' = -e cos p` on the branch.

Now write `p = pi - x_t`, `q = pi - x_u` (`x_t = t mod 2pi`,
`x_u = (t-beta) mod 2pi`), so `p, q in [-pi, pi]`,
`d := q - p in beta + 2pi Z` and `sin d = sin beta` on every branch.
Substituting the branch forms:

```
   Wtau = |p q| * ( sin d + d sinc(p) sinc(q) )        =: |pq| * fhat
   Wpot = |p q| * (-sin d + d sinc(p) sinc(q) )        =: |pq| * ghat
   Wwgt = |p q| *          d sinc(p) sinc(q)           =: |pq| * S
```

with `sinc(z) = sin z / z` **entire**.  The reduction has absorbed every
kink into a piece boundary: the four kink points are exactly `p = 0`
(`t = pi`), `q = 0` (`t = beta+pi`), `p = +-pi` (`t = 0`), `q = +-pi`
(`t = beta`) — the boundaries of the four sign quadrants of `(p, q)`.

Verified two ways in `validate.py`: symbolically with `sympy` on all **four
sign quadrants** (all twelve identities), and on **7144 interval samples**
where the enclosure of (folded mirror - reduced form) contains 0 with radius
at most `2.0e-57`.

### 2.3 The exact linear relation

Because `fhat + ghat = 2 S` identically,

```
    Wtau(beta,t) + Wpot(beta,t) = 2 Wwgt(beta,t)      for all beta, t
```

proved symbolically (`sympy`, exact) and verified in interval arithmetic
with residual radius at most **7.46e-58**.  The three determinants span only
a 2-dimensional space in this model too; any two determine the third, so
only `Wtau` and `Wpot` are analysed below.

### 2.4 The reflection symmetry

`h` is odd, `phi` and `|sin|` are even, hence `sA` is even, hence

```
   Wtau(2pi - beta, 2pi - t) = -Wtau(beta, t)      (same for Wpot, Wwgt)
```

and the mass direction `(-h(t-beta), h(t))` is negated, which is the
identity of the `(-,-)` fold — so **`y` is invariant**.  The whole boundary
set, and the whole census, is mirror-symmetric about `beta = pi`.  Measured
on the shipped classifier: identical at **1536 of 1536** grid teachers.  The
fundamental domain is therefore `beta in [0, pi]`.

---

## 3. The coincident stratum: root counts CLOSED, for every beta

The lattice structure first.  At `t = pi` **every atom of this kernel
vanishes**: `h(pi) = phi(pi) = |sin pi| = sA(pi) = 0`.  Consequently `t = pi`
and `t = beta + pi` are roots of all three determinants for every `beta`
(in the reduction, `p = 0` and `q = 0`, where the prefactor `|pq|`
vanishes), while `t = 0` and `t = beta` are roots of `Wwgt` only.  Interior
roots — roots with `|pq| != 0` — are exactly the interior zeros of `fhat`
(for `Wtau`) and `ghat` (for `Wpot`).

### 3.1 The analytic closure

`certify_analytic.py` (artifact `analytic_counts.json`) certifies the
interior-root counts for **every** `beta in (0, pi)` at once — no box
covering of the `beta` axis, hence no collars — by exact piecewise
inequalities on the reduced forms.  By the proved mirror symmetry (2.4) the
counts transfer to `(pi, 2pi)`.  The four smooth pieces of the `t`-period
are `L = (0, beta)`, `A1 = (beta, pi)`, `A2 = (pi, beta+pi)`,
`A3 = (beta+pi, 2pi)`; on `A1, A2, A3` the piece has `d = beta`, on `L` it
has `d = beta - 2pi`.

* **`Wtau` has no interior roots on `A1, A2, A3`.**  There
  `fhat = sin beta + beta * sinc(p) sinc(q)` with both sinc arguments in
  `[-pi, pi]`, where `sinc >= 0`; so `fhat >= sin beta > 0`.  Exact.
* **`Wpot` has no interior roots anywhere.**  On `L`,
  `ghat = -sin beta - (2pi - beta) S <= -sin beta < 0` (`S >= 0` there).  On
  `A1`/`A3` the sinc product is `sinc(x) sinc(x + beta) < sinc(beta)` for
  `x > 0` (sinc is strictly decreasing on `[0, pi]`), so
  `ghat < -sin beta + beta sinc(beta) = 0`.  On `A2` the arguments satisfy
  `a + b = beta`, and `sinc(a) sinc(b) > sinc(a+b)` — from the sympy-exact
  identity `(a+b) sin a sin b - ab sin(a+b) = a sin a (sin b - b cos b)
  + b sin b (sin a - a cos a)` and `sin x - x cos x > 0` on `(0, pi]` — so
  `ghat > 0`.  Hence **`Wpot`'s total count is 2 (the two lattice roots) for
  every `beta in (0, 2pi) \ {pi}`**: there is no potential curve in this
  model, the sharpest structural difference from the centered map.
* **`Wtau` on `L`: the count is decided at the piece midpoint.**  With
  `u = beta/2`, `v = pi - u`, `t = u + sigma`, the sinc product is
  `h(sigma) = (sin^2 u - sin^2 sigma)/(v^2 - sigma^2)`, even in `sigma`,
  and interior roots solve `F(sigma) := (2pi - beta) h(sigma) - sin beta
  = 0`.  The derivative numerator of `h` is governed by
  `Q(u, sigma) = sinc(2 sigma)(v^2 - sigma^2) - sin(u+sigma) sin(u-sigma)`,
  certified `> 0` on the whole triangle `0 <= sigma <= u <= pi/2` (15
  adaptive interval boxes, minimum certified margin `6.6e-2`, plus an
  analytic corner bound with margin constant `>= 0.55` covering
  `sigma > 1.171`), so `h` is strictly decreasing in `|sigma|`: `F` is
  maximal at the midpoint, where sympy-exactly
  `F(0) = (2 sin u / v)(sin u - v cos u)`, whose sign is the sign of
  `g(u) = sin u - (pi - u) cos u`; and `g` is strictly increasing on
  `[0, pi/2]` (`g' = 2 cos u + (pi-u) sin u > 0`) with `g(0) = -pi < 0` and
  unique zero `u* = beta*_1/2`.  At the piece endpoints `F -> -sin beta < 0`.
  Hence on `L`:

```
   beta < beta*_1 :  F < 0 everywhere            0 interior roots
   beta = beta*_1 :  F(0) = 0, F < 0 elsewhere   1 interior (double) root
   beta > beta*_1 :  F(0) > 0, one crossing
                     per side of the midpoint    2 interior roots
```

All seven algebraic identities behind this are sympy-exact
(`analytic_counts.json: checks.sympy`), including the identification of the
`L`-piece forms with the quadrant reduction of 2.2, and the one-dimensional
inequalities (`sin x - x cos x > 0`, the corner-bound constants, `g' > 0`)
are certified by series bounds plus single interval evaluations.  A float
sign scan at seven betas agrees with the analytic counts
(`checks.float_sanity`).

### 3.2 The counts

For `beta in [0, 2pi)`, with the exact degenerate cases split out, roots in
`t` over one `2pi` period are:

| | roots of `Wtau` | roots of `Wpot` |
|---|---|---|
| `beta = 0` | identically zero | identically zero |
| `(0, beta*_1)` | **2** (`t = pi`, `t = beta+pi`) | **2** (same two) |
| `beta = beta*_1` | **3** (a double root at `t = beta/2`) | 2 |
| `(beta*_1, pi)` | **4** | 2 |
| `beta = pi` | 2 (`t = 0`, `t = pi`); `sin beta = 0` makes all three determinants coincide, with zero set exactly the lattice | 2 |
| `(pi, beta*_2)` | **4** | 2 |
| `beta = beta*_2` | **3** (double root at `t = pi + beta/2`) | 2 |
| `(beta*_2, 2pi)` | **2** | 2 |

This closes what the first version of this certificate left open: the two
outer intervals `(0, beta*_1)` and `(beta*_2, 2pi)`, whose adaptive box
coverings had not terminated (the deciding margins degenerate like a power
of `beta` towards the column `beta = 0`, so a covering's cost diverges
there), and the six `beta` collars of total measure `4.0e-13` around
`beta*_1, pi, beta*_2`.  The box coverings that DID complete remain as
corroboration: 309 boxes on `(beta*_1, pi)` and 305 on `(pi, beta*_2)`, all
certifying `(Wtau, Wpot) = (4, 2)`, uncovered measure `2.0e-13` each
(`certificate.json: coincident_stratum`), and a 15-of-16-centre spot check
on the outer intervals, all `(2, 2)`, measure 0.517 = 11.6% of the two
intervals (`spot_check.json`; the 16th centre, `beta = 0.01`, sits too close
to the degenerate column for any box test).

### 3.3 `beta*`, exactly

The lens is born at the piece midpoint (`sigma = 0`, i.e. `t = beta/2`),
where `F(0) = 0` reduces to

```
   tan u = pi - u ,      u = beta/2 in (0, pi/2)
```

`g(u) = sin u - (pi-u) cos u` is strictly increasing (3.1), so the root is
unique; certified by interval bisection on the bracket `[1.11, 1.12]`:

```
   u*      in 1.11283481547935901488567225854
   beta*_1 = 2u*  in 2.22566963095871802977134451709
   beta*_2 = 2pi - beta*_1 in 4.05751567622086844715394224947
   beta*_1 + beta*_2 - 2pi  in [-1.044e-55, 1.044e-55]
```

This is the exact noncentered analogue of the centered
`tan u* = pi/2 - u*`.

### 3.4 The image of the coincident stratum in the map

| root | where | mass direction | `y` |
|---|---|---|---|
| `t = pi` | `beta` off the degenerate columns | `s1 = h(pi) = 0`, with nonzero `s0` | `y = 0` |
| `t = beta + pi` | `beta` off the degenerate columns | `s0 = -h(pi) = 0`, with nonzero `s1` | `y = +1` (= `-1`) |
| `t = 0` (`Wwgt` only) | `beta` off the degenerate columns | `s1 = h(0) = 0`, with nonzero `s0` | `y = 0` |
| `t = beta` (`Wwgt` only) | `beta` off the degenerate columns | `s0 = -h(0) = 0`, with nonzero `s1` | `y = +1` |
| lens pair | `beta in (beta*_1, pi) union (pi, beta*_2)` | generic | `y in (-1, 0)` |

Off the degenerate `beta` columns, the coincident image consists of **two
horizontal lines (`y = 0` and the `y = +-1` seam) and the torque lens**.
On the degenerate columns common-zero kernel vectors do not define `y`; those
columns are separate strata, with the displayed horizontal values only their
limits. The lens is
symmetric about `y = -1/2` under `t -> beta - t (mod 2pi)`, which sends the
mass ratio `rho` to `1/rho` and hence `y` to `-1-y`; both endpoints
`(beta*_1, -1/2)` and `(beta*_2, -1/2)` lie exactly on `y = -1/2`, and as
`beta -> pi` the two branches reach `y -> 0` and `y -> -1`, where they merge
with the lattice roots (at `beta = pi` the count drops from 4 to 2).  The
certified branch boxes are in `certificate.json:
coincident_stratum.beta_partition[*].boxes[*].lens_roots` (`y`-hull of the
tight boxes `[-0.5374, -3.3e-14]`, symmetry `y_upper + y_lower = -1` exact).

---

## 4. The separated stratum

### 4.1 The mass-free elimination

Two students at `th0 != th1` with weights `c0, c1`.  The two angular
(torque) equations solve the weights as `c0 = A(th1)/h(D)`,
`c1 = -A(th0)/h(D)` with `D = th0 - th1`, and the two radial equations are
the remaining balance

```
   F0 = kap c0 + phi(D) c1 - P(th0) ,   F1 = phi(D) c0 + kap c1 - P(th1).
```

`P` and `A` are linear in the teacher masses `(s0, s1)`, so a separated
critical configuration of a teacher with masses not both zero forces the
single mass-free equation `separatedDet = F0(e0) F1(e1) - F0(e1) F1(e0) = 0`
(`e0, e1` the unit teachers). Conversely, the centered theorem makes the
explicit first-row kernel vector satisfy both balances, but that vector may
be zero; it supplies a nonzero teacher only under a separate nonvanishing
check. A nonzero vector in the remaining rank-deficient cases requires a
separate kernel choice. **For the CENTERED kernel this elimination and
its converse are the Lean theorems `separatedDet_eq_zero_of_isSeparatedBalanced`
and `isSeparatedBalanced_of_separatedDet_eq_zero`
(`Planar/GeneralTeachers/SeparatedBoundary.lean`, stated with `kap = pi/2`,
`phiCos`, `couplingH`).  The noncentered instance used here — `kap = pi`,
`phiJ`, `HJ` — is the same two-line linear-algebra argument, re-derived in
this certificate and checked symbolically and numerically; it is not itself
a statement in the tree.**  The second-order analogue that IS in the tree
for this kernel is
`generalJKernelTeacherSchurDetT_eq_zero_of_separatedPair_degenerate`
(`Planar/SignedN2/FreeMassJ/NoncenteredGeneral/SchurBoundary/
SeparatedTypeBoundaryJ.lean`): the separated TYPE-change locus (vanishing
cleared Schur determinant) is mass-free.

The division-free form used for sweeping is `Dsep := h(D)^2 separatedDet =
a0 b1 - a1 b0` with the four column entries of `noncentered.sep_det`; `Dsep`
involves only `h` and `phi`, both `C^1`, so its value and gradient
enclosures are sound across the kink lattice.  Checked against the float
mirror of the division form on 3036 samples (max difference `4.06e-14`, the
mirror's own float error).  `Dsep` is antisymmetric under `th0 <-> th1` and
vanishes to fourth order on `D = 0`; the exclusion function is
`Psi = Dsep/h(D)^4`, which is `O(1)` as `D -> 0`.

There is a second clearing with a strictly better singular locus: solving
the RADIAL rows instead (`massDet = kap^2 - phi(D)^2 = pi^2 - phiJ(D)^2`,
which vanishes only at `D = 0 mod 2pi` — in particular not at the antipodal
gap `D = pi`, where `h` vanishes and the torque-row clearing degenerates)
gives `DsepWC`, and the sympy-exact factorisation

```
   DsepWC = massDet * Dred ,
   Dred   = hD^2 Wphi + hD (kap (Wp0 + Wp1) - phiD X) - massDet Wh
```

with five two-atom Wronskian groups (`Wp0`, `Wp1` are the potential
Wronskians at the two students; `Wphi`, `X`, `Wh` are cross-student
combinations; `noncentered.wc_columns`, `dred_grad_mv`; gradients validated
against central differences to `3.5e-39`).

### 4.2 The two-dimensional mass-free sweep: attempted, measured, not closed

The sweep runs interval branch-and-bound over the half domain
`(s, D) = (th1, th0 - th1)`, `s in [0, 2pi]`, `D in [delta, pi]` (the whole
stratum: of `D` and `2pi - D` exactly one is `<= pi`), excluding boxes whose
determinant enclosure omits zero and certifying graph boxes by interval
Newton.  It does not close, and the failure is now MEASURED from three
independent directions (thin `beta = 1.22`, `delta = 0.15` unless said):

* **Conditioning of the torque-row form** (`conditioning_measure.py` ->
  `conditioning.json`): uniform boxes of side `w = 0.1 / 0.05 / 0.025` give
  1827 / 7182 / 28728 boxes with 25.3% / 50.7% / 67.0% excluded by `Psi`.
  At the finest side the enclosure of `Dsep` is a median **3.7x** its true
  range on the box and the enclosure of `dDsep/dth0` a median **17.9x**
  (4104 sampled boxes, true ranges from 81 float samples per box), while
  only **3.2%** of cells genuinely contain the curve — the exclusion is an
  order of magnitude short of what the geometry allows.
* **The well-conditioned clearing does not help** (head-to-head, 400 random
  boxes per side): exclusion 22.2% / 52.2% / 69.0% (`DsepWC`) and 24.8% /
  50.5% / 66.2% (`Dred` with group-level mean-value refinement) against
  23.8% / 53.5% / 69.8% for the old form at `w = 0.1 / 0.05 / 0.025`.  The
  bottleneck is interval dependency in the four-column product, not the
  choice of clearing factor.
* **Newton never certifies the curve region**: at `w = 0.0125`, of 600
  random boxes 455 exclude, 0 pass either interval-Newton test (both
  partial-derivative enclosures straddle zero), 6 genuinely contain the
  curve, and 139 fail exclusion spuriously; the spurious failures do close
  under further bisection (3-255 extra evaluations each, measured on five),
  so the branch-and-bound converges — but along the curve nothing
  terminates it, and a full-domain run exhausts any budget tried (a
  400000-step budget at each of the 28 face witnesses returned the whole
  domain as undecided; `face_bb_attempt.json`).

**The honest verdict stands: the mass-free surface `{separatedDet = 0}` is
NOT certified as a complete enumeration**, in either clearing, at any
resolution reached. The per-teacher instrument of 4.3 replaces it: packed
branch-and-bound replays make the lists exact at twenty-two named teachers;
every other locate-then-certify list proves only existence and type of the
returned families (section 9.2).

### 4.3 The per-teacher instrument: locate-then-certify

For a FIXED teacher the separated solutions are isolated points, and the
conditioning is completely different.  Solve the weights from the RADIAL
rows, in the well-conditioned direction:

```
   massDet = kap^2 - phi(D)^2 = pi^2 - phiJ(D)^2
   c0 = (kap P(th0) - phi(D) P(th1)) / massDet
   c1 = (kap P(th1) - phi(D) P(th0)) / massDet
   G0 = c1 h(D) + A(th0) = 0 ,   G1 = c0 h(D) - A(th1) = 0
```

`phiJ` has range `[0, pi]` with `phiJ = pi` only at `D = 0`, so `massDet`
vanishes **only at student coincidence** — the antipodal stratum `D = pi` is
inside the certified domain, and separated families sitting exactly on it
are found (for example a saddle at `D = pi` for `beta = 2.6, y = -1/2`).
Until commit `2201dd0` the shipped `separatedScan` divided by `h(D)` with a
`|h(D)| > 5e-3` guard and was blind to that stratum; the shipped locator now
uses exactly this direction.

The enumeration is **locate-then-certify**: a dense float locator over
`(th1, D)` with `D` sampled GEOMETRICALLY from `1e-4` (a guess, never a
verdict), then an individual **2x2 Krawczyk existence-and-uniqueness test**
on a shrinking box around each candidate, in the `(s, D)` coordinates where
the cell is a genuine box.  Each certified family carries an enclosure
(typically `1e-12` to `1e-47` wide) and a certified `schurLabel` —
`widgets.js`'s own Schur test `T00 > 0 and det T > 0`, evaluated in interval
arithmetic; the interval entries equal the Lean
`generalJCramerSchurT00/T11/T01` divided by `massDet^2 > 0`
(`SeparatedTypeBoundaryJ.lean`), so the sign test is the tree's second-order
test (`GeneralJSeparatedPairCritical.lossSecondVariationNonneg_iff_schur`).
The unordered pair is de-duplicated modulo `2pi` and the EXACT FIT (both
students on the teacher lattice) is dropped, exactly as the classifier does,
because it is listed separately as `fit:global`.

**Two certificate KINDS, and what each does not give.**  This distinction is
the one a reader of the census map most needs, so it is stated here once and
every later claim keeps to it.

* An **existence certificate** is a Krawczyk box: `K(Z) subset int(Z)` on a
  box `Z` of student coordinates proves that a solution EXISTS in `Z` and is
  UNIQUE there.  It is a statement about `Z` alone.
* An **exclusion certificate** is a signed residual enclosure on a box: if
  the outward-rounded enclosure of `G0` (or of `G1`) over the box does not
  contain zero, the box holds NO solution.  Exclusion over a region is what
  proves nonexistence there.

**A pile of existence certificates is not a completeness proof.** Every
separated family returned by the locator is individually certified to exist,
be locally unique, and have the stated Schur type; that fact alone does not
exclude an unseeded family. Completeness requires an exclusion cover over the
whole complement of the certified boxes. The packed branch-and-bound replays
now provide such a cover at twenty-two named teachers, totaling 2,368,457
leaves with no failed or undecided leaf. At those teachers the lists are exact.
At every other teacher discussed through locate-then-certify alone, the honest
reading remains “at least the listed families.” Pointwise completeness does not
by itself propagate to a proposed cell or to an interpolated wall arrangement.

**The geometric `D`-grid is not a detail.**  The separated families that
carry the census near `beta = pi` sit at student gaps `D ~ 0.02 - 0.07`, and
those near the strata `y = 0`, `|y| = 1` at gaps proportional to the
distance from the stratum (section 6.3).  A uniform `D`-grid misses them
entirely — that was the first version of this scan, and, until `2201dd0`,
the shipped map's locator.

---

## 5. The census IS the sign chart, plus a separated term

At a torque root where `(-h(t-beta), h(t))` is nonzero, the mass direction is
`(s0,s1) = kappa (-h(t-beta), h(t))` for a real `kappa != 0`, hence

```
   P(t)   = s0 phi(t)  + s1 phi(t-beta)  = -kappa Wpot(beta,t)
   W(t)   = s0|sin t|  + s1|sin(t-beta)| = -kappa Wwgt(beta,t)
   tau(t) = P(t) - 2 W(t)                = +kappa Wtau(beta,t)
```

the last using `Wtau + Wpot = 2 Wwgt` (2.3).  `kappa^2` cancels in
`widgets.js coincidenceTypeAt`, so

```
   a trap is present at this root   <=>   Wtau * Wpot < 0
   its label is  @positive          <=>   Wtau * Wwgt > 0     (else @mixed)
```

verbatim the centered identity — the only input is `K'' + K = 2|sin|`,
which both kernels obey.

**But the noncentered census is NOT a function of the three determinants
alone.** Unlike the centered classification, the noncentered one enumerates
separated families and emits `separate:trap@positive` /
`separate:trap@mixed` rows. Every pointwise row list below therefore has
two halves: the sign chart of section 3 (complete torque-root isolation plus
certified determinant signs), and the separated family list of section 4.3.
The latter is exact at the twenty-two replayed exhaustive teachers and is only
a list of certified existing families at every other locate-then-certify point.

---

## 6. The 1-dimensional and 0-dimensional pieces

Map coordinates `x = beta in [0, 2pi)`, `y in (-1, 1]`.

### 6.1 One-dimensional

**1D-1 / 1D-2. torque lens, upper and lower branch.**  `Wtau = 0` interior
roots, `beta in [beta*_1, beta*_2]`; the branches are mirror images under
`t -> beta - t` (`y -> -1-y`), meet at the two double roots
`(beta*_i, -1/2)`, and reach `y -> 0` / `y -> -1` as `beta -> pi`.

**1D-3 .. 1D-6. the lines `y = 0` and `y = +-1`.**  Images of the lattice
roots `t = pi, beta+pi` (all three determinants) and `t = 0, beta` (`Wwgt`
only); equally the degenerate strata `s1 = 0` and `s0 = 0`.

**1D-7. degenerate column `beta = 0`.**  All three determinants vanish
identically; coincident teacher directions.

**1D-8. degenerate column `beta = pi` (antipodal).**  `sin beta = 0`, so all
three determinants coincide, with common zero set the lattice `{0, pi}`;
here the lens branches merge into the lattice roots.

**1D-9. the `h(D) = 0` locus of the separated stratum** (`D = 0 (mod pi)`),
a locus in `(beta, th0, th1)`, not a curve of the map; `D = 0` is the
modelling boundary of the separated stratum, `D = pi` is interior to the
certified domain of 4.3 (only the pre-`2201dd0` classifier excluded it).

**Proposed traces 10 / 11: sampled candidate transitions in the same-sign sector.**
Across the short `y`-windows returned by `ywall`, the
locate-then-certify family lists differ. Every returned family is individually
certified, but these teachers have no exclusion cover for the separated
complement. The samples therefore do not certify a count drop, an
annihilation, or a census boundary. The renderer joins them into a proposed
trace `w_a(beta)`; trace 11 is its pointwise image under
`y -> -1-y`. No continuous global wall arc is enclosed here.

**Proposed traces 12 / 13: sampled candidate transitions in the mixed sector.**
These organize the additional collar row of 6.3; they are not certified
one-dimensional census pieces.

### 6.2 Zero-dimensional

The certified zero-dimensional pieces are the two lens double roots
`(beta*_1, -1/2)`, `(beta*_2, -1/2)` (exact
`y = -1/2`: `h(t) = -h(t-beta)` there); the two lens terminations on
`beta = pi` (`y = 0` and the seam, genuine `0/0` mass directions); the
corners of the degenerate columns with `y = 0` and the seam. The
four sampled candidate-trace/column meetings are vertices only in the proposed
incidence graph of section 7. The lens stays strictly inside
`y in (-1, 0)` for `beta != pi` (at a lens root `h(t)` and `h(t-beta)` have
opposite signs), so it meets the horizontal strata only at its `beta = pi`
terminations.

### 6.3 The mixed collar: an additional observed selected-row list

At representative points in the sector `y > 0`, the proposed-label row list
of section 7 carries two
`separate:trap@mixed` rows. In collar samples hugging the strata
`y = 0` and `|y| = 1`, the locator returns one:

> **Locate-then-certify witnesses.** At
> `(beta, y) = (0.15, 0.004)` and the teacher-swap mirror
> `(0.15, 0.996)`, the returned list is
> `coincident:trap@positive x2 | fit:global | separate:trap@mixed`.
> The coincident sign chart is complete and every located separated family is
> Krawczyk-certified (`wall_pos.json`), but no exclusion cover disposes
> of unseeded separated families. The committed `map_scan` supplies the
> same evidence class at `(0.05, 0.02)` and `(0.05, 0.98)`.

The proposed mechanism is represented by candidate traces 12/13: the second located family's
student gap `D_2` shrinks with `y` (measured
`D_2 ~ 0.5 y` at `beta = 0.75`). The pointwise first-success evidence is
stored in `wall_pos.json` and re-audited in `jacobian_fix_delta.json`.

The table below records the **post-Jacobian-fix** pointwise scope (commit
`5e5adf1`, 9.5). “Certified” modifies each returned family, not
completeness of the list. Near the proposed transition the corrector usually
declines on the one-trap side, so only the first row is a two-sided sampled
pair; the other rows are one-sided first-success heights.

| `beta` | certified located one-trap witness | first certified located two-trap height | scope of the sample |
|---|---|---|---|
| 0.15 | `y = 0.004` | `y = 0.200` | two endpoint witnesses; no intervening claim |
| 0.30 | -- | `y = 0.080` | one-sided first success only |
| 0.50 | -- | `y = 0.020` | one-sided first success only |
| 0.75 | -- | `y = 0.012` | one-sided first success only |
| 1.00 | -- | `y = 0.008` | one-sided first success only |
| 1.50 | -- | `y = 0.005` | one-sided first success only |
| 2.20 | -- | `y = 0.005` | one-sided first success only |
| 2.60 | -- | `y = 0.005` | one-sided first success only |
| 3.00 | -- | `y = 0.020` | one-sided first success only |

Both moves are in the same direction — with the corrected Jacobian the SECOND
located family certifies at smaller `y`, so the sampled first-success
heights are lower than before. They are nonmonotone across the sampled gaps,
falling through the middle values and rising again at `beta = 3.00`.
This does not prove that one wall curve exists between samples or is itself
nonmonotone.

The rerun changes two sampled first-success heights. It does not certify a continuous
stratum, a `beta`-extent, a one-dimensional piece count, or a nine-region
global refinement. By the exact `y -> 1-y` symmetry the sampled evidence has
mirror points below `|y|=1`; at `(0.15,0.004)` and its mirror the
listed separated row is certified to occur, but the list is not exhaustive.
Whether boundary curves exist, how any such curves continue between the nine
sampled gaps, and the collar interior remain open; most near-transition points on the
one-trap side make Krawczyk decline.

This finding is POST-locator-fix: the pre-`2201dd0` classifier could not
see the collar at all (its uniform seed lattice missed every near-coincident
family), and the first version of this certificate saw it only as an
unresolved anomaly at `beta = 0.05` (then section 9.5).

---

## 7. The proposed arrangement and its combinatorial count

**Domain.**  `beta in [0, pi]` is a fundamental domain (2.4); `y` is a
circle of circumference 2; so the domain is a closed cylinder whose boundary
circles are the columns `beta = 0` and `beta = pi`, both genuine strata.
A second exact symmetry — swapping the two teacher atoms, a rotation
composed with the `beta`-mirror — acts at fixed `beta` by `y -> -1-y`,
making the map symmetric about `y = -1/2` in the same-sign sector and about
`y = +1/2` in the mixed sector (the pairing F2/F6, F3/F7 below, and the
two collars of 6.3).

**Wall evidence.** The torque-lens root count, the horizontal strata
`y = 0` and `y = +-1`, and the two boundary columns have the analytic or
interval scope stated in sections 2--3. The same-sign-sector traces 10/11
come from individually certified returned families at the sampled windows of
section 6; the lists are not exhaustive. Joining those samples into global
arcs is an interpolation. If those proposed arcs
have the displayed incidences and no additional wall occurs, cutting the
cylinder gives

```
   vertices  V = 11 ,  edges  E = 18 ,  faces  F = 7
   V - E + F = 0 = chi(annulus)                    [conditional Euler check]
```

The proposed eleven vertices are the lens double root, the two lens
terminations on
`beta = pi`, the four column/line corners, and the four wall/column
meetings; the eighteen edges: 2+2 lens branches each cut once by a proposed
candidate-transition arc, 1 + 1 horizontal lines, 2+2 proposed transition
arcs each cut once by the lens, and
4+4 column arcs.  Under the same assumptions the full-torus count would be
`V = 14, E = 28, F = 14`.  These are combinatorial bookkeeping values, not
a certificate that the proposed graph is the complete global arrangement.
The mixed collar traces of 6.3 are not in this count: one gap has a sampled
transition pair and eight have one-sided first-success heights, so the
suggested refinement to **9** cylinder regions (**18** on the torus) is
likewise not certified globally.

**Seven proposed labels and their representative points.**  The labels
F1--F7 name regions in that proposed graph; they do not assert that a global
cell has been enclosed.  All 28 entries of `faces_fixed.json` have
locate-then-certify existence boxes and agree with the shipped classifier at
commit `2201dd0`.  Separately, the packed exhaustive branch-and-bound replays
of section 9.2 close twenty-two named teachers, including every representative
in the table below.  Thus the displayed census is exact at the displayed
point, but pointwise completeness does not prove constancy over the proposed
region.  The Jacobian repair of section 9.5 moved the F3 and F7 representatives
from `(3.13, -0.08)` and `(3.13, -0.92)` to `(3.08, -0.08)` and
`(3.08, -0.92)`.

| proposed label | proposed region | representative point | `n` | exact census at that replayed teacher |
|---|---|---|---|---|
|  F1 | mixed sector `0 < y < 1`, outside the collars of 6.3 | `(0.75, 0.40)` | 2 | `coincident:trap@positive \| coincident:trap@positive \| fit:global \| separate:trap@mixed \| separate:trap@mixed` |
|  F2 | `w_a(beta) < y < 0`, outside the lens | `(0.75, -0.02)` | 2 | `coincident:trap@mixed \| fit:global \| separate:trap@positive` |
|  F3 | `w_a(beta) < y < 0`, inside the lens | `(3.08, -0.08)` | 4 | `coincident:trap@mixed \| coincident:trap@mixed \| fit:global \| separate:trap@positive` |
|  F4 | `-1-w_a < y < w_a`, outside the lens | `(0.75, -0.40)` | 2 | `coincident:trap@mixed \| fit:global` |
|  F5 | `-1-w_a < y < w_a`, inside the lens | `(2.60, -0.50)` | 4 | `coincident:trap@mixed \| coincident:trap@mixed \| fit:global` |
|  F6 | `-1 < y < -1-w_a`, outside the lens | `(0.75, -0.98)` | 2 | `coincident:trap@mixed \| fit:global \| separate:trap@positive` |
|  F7 | `-1 < y < -1-w_a`, inside the lens | `(3.08, -0.92)` | 4 | `coincident:trap@mixed \| coincident:trap@mixed \| fit:global \| separate:trap@positive` |

`w_a(beta)` is proposed trace 10 (sampled locate-then-certify windows
in 6.1), and
the sampled brackets place the proposed crossing `beta_c` in
`(2.60, 3.05)`: at `beta = 2.60` the certified lens span is
`y in (-0.72, -0.28)`, while `(3.05, -0.09)` has individually
certified returned rows for the overlap label. This supplies nonempty witness
sets; it does not enclose their full `beta`-extent or certify the interpolated
crossing arc.

**Observed pointwise census values.** The seven representatives, all among the
twenty-two exhaustive replays, display five exact values (F2/F6 and F3/F7
agree by the `y -> -1-y` symmetry). The collar witnesses of 6.3 supply
an additional row pattern whose listed separated family is certified to
exist, not a sixth exact census. Calling these the exhaustive set of global
region values would require the missing arrangement proof.

**The shipped map.**  At commit `2201dd0` the noncentered map draws 49
components at depth 8 (19 below the naming threshold, 0.0034% of the
diagram); the big components correspond visually to the proposed labels
doubled by the
`beta`-mirror, the collar strata of 6.3 at small `beta`, and the rendered
measure-zero lines.  At commit `f430fab` it drew 1898 (1868 sub-threshold,
speckled across the antipodal band) -- section 8.  These component counts are
measurements of a renderer/classifier, not a proof of the arrangement.

---

## 8. The antipodal band: the defect, and the fix, against named commits

**The measurement that identified the defect** (against `widgets.js` at
commit `f430fab`; artifact `band_compare.json`).  On a `12 x 12` grid over
`beta in [3.05, 3.235]`, `y in (-1, 1)`, plus a partial refinement over
`beta in [3.128, 3.156]` — 158 points — the then-returned row lists
(sign chart plus Krawczyk-certified located separated families) were compared
with the shipped
`censusSignature`, loaded through the generator's own DOM shim:

```
   points 158    certified 146    shipped agrees 94 of the 146
   disagreements 52, EVERY one with the shipped classifier reporting FEWER
   separated traps than are certified to exist (plus 5 more at points the
   certificate declines to certify)
```

The returned lists passing the then-current pointwise tests were constant
across the sampled `beta` values — three lists depending on `y` alone
(the same lists as at `beta = 0.75`) — while the
shipped signature took five values in a speckled pattern; that speckle is
what rendered as 1868 fragments.  A one-sided error at 52 of 52 certified
disagreements is a bug signature, not sampling noise, and three mechanisms
were identified, all in the locator, none in the tolerance: the
`110 x 110` uniform seed grid stepped over the near-coincident families
(gaps `D ~ 0.02-0.07` — one grid cell) that carry the band's census; the
`|h(D)| > 5e-3` guard blinded it to the antipodal stratum `D = pi`; and
`findRoots` missed one of four torque roots within `~3e-3` of `beta = pi`.

**The fix** (commit `2201dd0`, `separatedScanJ`): the shipped locator now
uses this certificate's well-conditioned direction (4.3) with a geometric
gap grid and W-graded seeding, no absolute pre-filter, and acceptance
against the true balance.  Re-measured on the same 158 points
(`band_comparefix.json`): **146 of 146 certified points agree**; the 10
remaining disagreements all sat at the 12 points the certificate then
declined to certify (`|beta - pi| < 0.0155`), where the shipped classifier
reports MORE separated families than the certificate's partial list — the
expected direction at points where the certificate's own locator is below
its resolution.  The shipped band census takes 4 values (was 6), and the
map's noncentered block dropped from 1898 components to 49.

**The Jacobian fix supersedes those numbers** (commit `5e5adf1`, 2026-08-08;
see 9.5).  Until then `sep_res` computed one block of the Krawczyk Jacobian
with a wrong sign, so the interval `dG0/d(th1)` was an INVALID enclosure and
the Krawczyk verdicts on this band were unreliable in both directions.  The
band was re-certified with the corrected Jacobian at the scan's own
`N = 180` (`band_rescan_fixed.py` -> `bandfix2_*.json`,
`band_comparefix2.json`):

```
   points 158   evaluated 156 (2 time out, see below)
   certified 137    shipped agrees 137 of the 137
   disagreements 19, EVERY one at a point the certificate DECLINES
```

The structural claim survives verbatim — **every disagreement still sits at
a point the certificate itself declines to certify** — but the headline is
`137/137`, not `146/146`, and the declined set grows from 12 points to 19.
All 19 lie in the antipodal pocket `|beta - pi| in [0.0068, 0.0240]`, which
REACHES OUTSIDE the `< 0.0155` zone this section previously declared; the
scope that the grid supports is stated in 9.5.

The two points that time out are `beta = 3.129` and `3.131` at `y = 0`, and
they leave **no hole in the certified table**: `y = 0` is the mass stratum
itself, where `masses_at` returns `s1 in [-6.3e-49, 7.4e-49]`, i.e. `s1 = 0`
and a ONE-teacher configuration outside the two-teacher model this
certificate is about.  Both were already `certified = False` BEFORE the fix
(0 certified torque roots, 91-100 undecided torque boxes) with the
degenerate `n_separated = 4996` — the signature of `s1 -> 0`, where the
student load collapses to `s0 * num(D, .)` and its zero set becomes the
whole teacher lattice — and neither appeared among the pre-fix
disagreements.  They are therefore excluded from both the `146` and the
`137`, and the stall is in the separated half only (`torque_roots` returns
in 2.6 s).  Their presence in the 158-point grid is an artifact of the grid
having landed exactly on a 1-D stratum.

The `SEL_REL` tolerance change that fixed the centered map's slivers
correctly changed nothing here: this was a defect of the locator, not of the
tolerance.

---

## 9. The unresolved set

### 9.1 Coincident stratum: nothing

The analytic closure of section 3 leaves no collar and no uncovered
interval on the coincident stratum.  The columns `beta = 0` and `beta = pi`
are exact degenerate cases, not gaps (all three determinants vanish
identically at `0`; common zero set = lattice at `pi`).

### 9.2 Completeness of the separated locator

The packed exhaustive branch-and-bound replays now cover twenty-two named
teachers across F1--F7.  Together they contain 2,368,457 leaves, with zero
failed leaves, zero undecided boxes, and certified disjointness of the
enclosed classes.  At those twenty-two teachers the separated family lists
are exact.  The table and replay commands in the companion README are the
authoritative ledger.

Outside those named teachers, a locate-then-certify census still begins with
a floating-point locator (`N = 160--260` angular seeds and about 240 geometric
gap values from `1e-4`).  Every returned family is then certified
individually, but an unseeded family is not excluded.  The exact point
enumerations do not establish a count over any proposed F1--F7 region; that
regional transport requires the boundary-free and fold-free hypotheses of
the Lean count-constancy schema.  The mass-free sweep discussed below does
not discharge those hypotheses globally.

**How the pointwise closure was reached, and what it does not establish.**
The strict
SIGN LAW of the separated scalar is now a Lean theorem
(`Planar/SignedN2/FreeMassJ/SeparatedCount/SignLawJ.lean`), together with
its geometric consequence over every teacher with `s0*s1 != 0`: each
student's two teacher-offsets STRADDLE the gap lattice for positive masses
(`LatticeStraddleJ.lean`) and lie on the SAME open side for mixed signs
(`MixedSameSideJ.lean`).  That dichotomy is a box-exclusion test of four
endpoint comparisons — no interval evaluation of the residual, hence none of
the dependency blow-up that stalled the sweep of 4.2 — and it prunes
`50-91%` of the `(th1, D)` torus for positive teachers, `5-49%` for mixed
ones.  With it the per-teacher branch-and-bound reaches a residual undecided
area of `3.1e-3` (F4 witness, 400000 steps) against `18.68` — 95% of the
domain — for the same budget without it (`face_bb_signlaw.py`).
**It now TERMINATES, and at the F5 representative it CLOSES.**  At
`(beta, y) = (2.6, -1/2)` the search ends with

```
   complete = true    undecided boxes = 0    undecided area = 0
```

— the entire `(th1, D)` domain is either excluded by a certificate or
Krawczyk-certified, with no residual.  The certified dedup returns 5 raw
boxes, one certified merge, 4 classes with disjointness certified, of which
one is the EXACT FIT and **three are separated families**, at
`D = [1.8836, 1.8836, pi]` — exactly the census's
`[3.141593, 1.883569, 1.883569]`.  So at that witness the census carries
EXACTLY three separated families, not merely at least three: **the residual
obligation of this section is discharged at that point.**  It is a POINT and
not a face: complete enumerations at teachers carrying the proposed F5 label
have different counts. Thus that proposed label cannot be one count region;
the drawing's subdivision in `f5_count_band.json` is prospecting, not a
certified face subdivision.

Four defects had to be repaired, and three were the SAME defect — *a solution
lying exactly on a boundary of the decomposition can never be
Krawczyk-certified, because the test needs it strictly interior, while the
cells touching it can never be excluded either, so they subdivide forever.*
The exact fit sits at the seam `th1 = 0`; the antipodal family sits at
`D = pi` exactly; and the exact fit's gap sits at `D = beta`, an atom line.
The repairs are a generic seam offset, `D` extended past `pi` through the
`D <-> 2pi - D` fold (which over-counts rather than misses, and whose extra
copies the dedup removes by that same involution), and each atom-derived cut
replaced by a PAIR straddling it, so the atom line lies strictly inside a thin
cell.  The fourth defect was the ORDERING: best-first by area refines the whole
domain to certification scale before reaching any solution, so it certifies
nothing; depth-first pays that cost only where it is needed.

The lever that made it cheap is the MEAN-VALUE form
`G(X) subset G(xc) + J(X).(X - xc)`.  Plain interval evaluation of the
residual loses a factor 3.7 to dependency (4.2), so on small boxes it cannot
sign a residual that is nowhere near zero.  Measured at the same witness, the
centered form performs **24839 of the 24972 exclusions** and cuts the run from
464706 steps with 1663 residual boxes to 53172 steps with 2.  It is available
only because the Jacobian is now correct — the sign error of 9.5 is exactly
what made this form unusable before.

This first F5 closure was subsequently extended to all twenty-two teachers
listed in the README.  It remains pointwise: no global F1--F7 region count
follows, and the older `count_match.json` face fields are not a current
completeness ledger.

### 9.3 The mixed collar samples

The additional collar row is certified to occur at the witnesses of 6.3, but
the separated lists there are not exhaustive. Only `beta = 0.15` has a
sampled transition pair; the other eight gaps have one-sided first-success
heights. The scan reaches `beta = 3.00`, but these observations do not
locate the functions `w^+(beta)` and `w^-(beta)` or prove global
continuation. The near-transition strips where Krawczyk declines remain outside
that pointwise evidence.  The same-sign scan likewise has an uncertified
strip near `y = 0` at small `beta` (`certificate.json: map_scan`), while the
antipodal-band re-audit and its precise scope are recorded in 9.5.

### 9.4 The per-teacher coincidence collar

The per-teacher enumeration certifies families at gaps `D >= dmin = 1e-4`
and treats `D < dmin` as coincident territory; a separated family with
`D < 1e-4` would be missed.  The collar families found near the strata
(6.3) show `D` scaling linearly with the distance to the stratum, so this
collar is the `y ~ 2e-4` neighbourhood of the strata — far below every
witness used.

### 9.5 The `sep_res` Jacobian sign error (found and fixed 2026-08-08)

Until commit `5e5adf1` the Krawczyk Jacobian in `census_cert.sep_res`
carried

```
dc1_1 = (kap*A1 + hD*P0)/M + c1*q/M          # WRONG
```

where `D = th0 - th1` and `phi' = -h` give `d(phiD)/d(th1) = +hD`, hence
`dc1_1 = (kap*A1 - hD*P0)/M + c1*q/M`.  The consequence was not a wide
enclosure but an **invalid** one: the interval `j01 = dG0/d(th1)` failed to
contain the true derivative range in **188 of 300** random boxes spanning
`beta in (0.2, 3.0)`, `y in (-0.9, 0.9)`.  An invalid derivative enclosure
breaks the Krawczyk test's hypothesis, so every `unique` / `empty` verdict
computed with it is void.  After the fix: 0 invalid enclosures in 300 boxes,
all four entries (`sep_jacobian_audit.py`, the permanent regression guard).

**What it did NOT touch.**  The residual itself is correct — `sep_res` and
the independent float mirror in `count_scan.py` agree to 9 digits — so the
float locator, which never forms this Jacobian, is unaffected; so is the
whole coincident half (section 3, analytic closure), and so is every Lean
theorem, none of which consumes it.

**What it did touch, measured item by item** (`jacobian_fix_delta.py`,
`jacobian_fix_delta.json`; 301 stored points re-run, 26 changed):

The JSON retains `collar_wall` and `fold_wall` as legacy set keys for artifact
compatibility. Its `legacy_set_names` field gives the current pointwise/sampled
meaning; neither key asserts a certified wall.

| set | points | changed |
|---|---|---|
| face witnesses (7.x) | 28 | 2 |
| mixed-collar witness scan (6.3) | 45 | 3 |
| mixed-sector candidate-transition scan (6.3) | 70 | 2 |
| antipodal band (8) | 158 | 19 (17 evaluated + 2 TIMEOUT) |

Every change lies in one of exactly two pockets, and the bulk of the stored
witness set is unchanged:

* **Pocket A, the antipodal pocket.**  All 19 band changes and both face
  changes.  Post-fix, every DECLINED band point has `|beta - pi| <= 0.0240`;
  the next distance the grid samples is `0.0376`, and all 84 points at
  `|beta - pi| >= 0.0376` are certified AND agree with the shipped
  classifier.  So the scope this grid supports is **`|beta - pi| >= 0.0376`
  certified**, against the `< 0.0155` no-go zone section 8 previously
  declared — the old zone was too narrow by at least `0.0085`.  Inside the
  pocket certification is mixed, not uniformly absent.
* **Pocket B, the `y > 0` collar strip**, `y ~ 0.008-0.02` at
  `beta in {0.5, 1.0}`: 3 collar and 2 mixed-sector
  candidate-transition-scan points.  Here the fix ADDS
  a separated family and certification, lowering the sampled first-success
  heights (table in 6.3).

The two directions matter: near `beta = pi` the bug FABRICATED families and
certification, on the collar strip it SUPPRESSED them.  It was not a
one-sided over-certification, so no "the old numbers were conservative"
reading is available.

**The two invalidated face witnesses were relocated, not deleted**
(`faces_fixed.json`, successor to `faces_0.json`, provenance kept): F3's
`(3.13, -0.08)` and F7's `(3.13, -0.92)`, both at `|beta - pi| = 0.0116`,
are marked `superseded_by_jacobian_fix` with their pre-fix records, and
replaced by `(3.08, -0.08)` and `(3.08, -0.92)` at `|beta - pi| = 0.0616`.
Each replacement reproduces the census string and Schur pattern assigned to
its proposed label. The table again carries 28 pointwise-certified witnesses,
with one observed census string per label; the minimum
`|beta - pi|` over all 28 witnesses is `0.0416`. This relocation repairs the
point table and makes no global-cell claim.

**How it was found, and why nothing else caught it.**  A certified
hull-merge dedup asked Krawczyk to re-verify its own stored enclosures and
got `empty` on boxes the same code had just certified `unique` — two
verdicts that cannot both hold.  No existing check could have seen it:
residual-level validation passes (the two implementations agree), and
precision sweeps say nothing, because outward rounding only makes Krawczyk
FAIL, never succeed spuriously.  Only an enclosure that is WRONG can
manufacture a verdict, and only a check of the interval DERIVATIVE against
sampled finite differences of an independent implementation can see that.
That check is now `sep_jacobian_audit.py` and should be run after any edit
to `sep_res`.

### 9.6 Why a box enumeration fails to close: two defect classes

Closing the per-teacher enumeration took five repairs.  Four of them fall into
two classes that recur, and both are worth stating separately because they
present identically — the search runs to its budget with a small residual —
and have completely different fixes.

**Class A: a solution lies exactly on a boundary of the decomposition.**  The
Krawczyk test certifies a solution only when it is STRICTLY interior to the
box; a solution on a cell boundary is interior to neither neighbour, and the
cells touching it can never be excluded either, because the residual genuinely
vanishes there.  So they bisect forever, and the residual is a pair of adjacent
cells sharing the offending boundary.  Three instances here, each found only
after the previous was fixed:

1. the EXACT FIT at the seam `th1 = 0`, the edge of the `s`-period;
2. the ANTIPODAL family at `D = pi` exactly, the edge of the `D`-range — and
   then, after the range was padded, a second time because `pi` was still
   being inserted as a `D`-CUT;
3. the exact fit's gap at `D = beta`, an ATOM LINE, which the cut set puts on
   a boundary by construction.

The fixes are all the same shape: move the boundary off the solution.  A
generic seam offset; `D` extended past `pi` through the fold `D <-> 2pi - D`
(which duplicates near-antipodal families rather than missing them — the safe
direction — with the duplicates removed afterwards by that same involution);
and every atom-derived cut replaced by a PAIR straddling it, so the atom line
lies strictly inside a thin cell.  **A boundary that carries mathematics is
the wrong place to cut.**

**Class B: a threshold set above the scale the mathematics needs.**  Distinct
from class A, and it was misdiagnosed as class A at first because the residual
looks the same — two adjacent minwidth cells bracketing a known solution.  The
cause is that Krawczyk certifies only inside a CONTRACTION RADIUS fixed by the
conditioning of the Jacobian at that solution.  At the F2 family
`D = 0.35938` the `(s, D)` Jacobian determinant is `-5.5e-3`, and the test
returns `unique` at radius `1e-5`, `1e-6`, `1e-7` but `None` at `1e-4` and
above.  The minimum box width was `2e-5`, ABOVE that radius, so the search
stopped bisecting before the test could ever fire and the family was silently
lost — the certified count came out one short of the census at exactly the
witnesses carrying an ill-conditioned family.

> **The minimum box width must sit below the contraction radius of the
> worst-conditioned solution in the domain, or that solution is lost.**  A
> minwidth residual bracketing a known solution is the signature, and the
> check is one line: run the test on shrinking boxes about the known solution
> and read off where it starts succeeding.

Lowering `minw` from `2e-5` to `1e-8` closed the witness with no other change
(`certified_via_inflation = 0`, i.e. the ordinary test fired once the boxes
were small enough).  The same mis-scaling had defeated the epsilon-inflation
fallback added for class A, whose ladder started at `1e-4` — above the
contraction radius, so it could never fire; it now starts at `1e-9`.

The fifth repair was neither: it was the ORDERING (9.2), where best-first by
area refines the whole domain to certification scale before reaching any
solution.  Depth-first pays that cost only where it is needed.

### 9.7 The completed 68-tile regional atlas

One regional statement is now closed without asserting the proposed global
arrangement.  On the strip
`censusStripJ 0.137 0.001 3.13`, the closed teacher rectangle
`beta in [0.5,2.2]`, `y in [-0.45,-0.1]` is covered by 68 tiles.  Their 68
merged regularity forests total 41,371,573 leaves and 206,853,513 preorder
nodes. The 272 physical-boundary forests add 41,724 leaves and 66,040 nodes,
for all-forest totals of 41,413,297 leaves and 206,919,553 nodes, with zero
undecided leaves and zero unresolved volume.  The witness has full zero
count `N = 2` and selected-minimum count `N = 1`.  The manifest is
`minimum_atlas_f4_census_safe.json`; the detailed replay ledger is
`MINIMUM_DFS.md`.

The Lean declarations
`noncentered_census_map_zeroCount_constant_on_atlas` and
`noncentered_census_map_selectedMinimumCount_constant_on_atlas` prove the
propagation result for an arbitrary atlas satisfying their hypotheses.  They
do not contain the 68 concrete tiles or prove the replay obligations.  Those
obligations are checked externally, and the identification of the checker's
interval expressions with the Lean constants `censusAngleMapJ`,
`censusCramerSchurT00J`, and `censusCramerSchurDetTJ` is human-checked rather
than a Lean theorem.  The regional conclusion therefore rests on the external
replay plus that stated correspondence, followed by the conditional Lean
theorem.  It says nothing about the global F1--F7 arrangement outside this
rectangle.

---

## 10. Honesty notes

* Every count in section 3 holds for **every `beta`**, by the analytic
  closure; the box coverings and the spot check are corroboration, not the
  evidence.  The closure's own instruments are exact sympy identities plus
  interval-certified one- and two-dimensional inequalities
  (`analytic_counts.json`).
* The coincident classification of section 3 is certified-complete.  For the
  separated stratum, every family listed by locate-then-certify is proved to
  exist, be locally unique, and have the stated Schur type.  At the twenty-two
  named teachers in the replay ledger, the exhaustive exclusion cover also
  proves that the list is complete.  At other isolated teachers it still means
  “at least the listed families.”  Neither kind of point evidence proves a
  global F1--F7 cell census.
* The 68-tile rectangle of 9.7 is a distinct regional certificate: external
  interval replay discharges the hypotheses of a generic Lean propagation
  theorem, subject to the explicitly human-checked checker-to-Lean expression
  correspondence.  The generic theorem itself does not encode the concrete
  tiles or their counts.
* The claim that the three determinants exhaust the census boundary is
  **not** made here; the separated stratum has sampled candidate-transition
  traces, but the
  locate-then-certify windows do not certify them as census boundaries. The mass-free
  elimination used for the separated sweep is Lean-proved for the CENTERED
  kernel and re-derived here for this one (4.1); the Schur type test equals
  the tree's second-order test up to a positive factor (4.3).
* Validity domain: `beta in [0, 2pi)` with fundamental domain `[0, pi]`
  (the symmetry of 2.4 is proved from the kernel parities), `t` in one
  `2pi` period, noncentered kernel only.  All widget comparisons name their
  commit: `f430fab` (pre-fix) or `2201dd0` (fix).
* The shipped mosaic is a rendering of `widgets.js` at its finite resolution,
  at both commits.  What section 8 reports broken-then-fixed is the
  classifier's locator, not the generator or renderer.  Pixel components are
  empirical output and do not certify the topology of the continuous wall
  arrangement.

---

## Files

```
  noncentered.py             kernel module: float mirror of widgets.js + mpmath.iv
                             interval forms, both mass-free clearings
  ref_widgetsJ.js            reference generator (literal widgets.js functions,
                             identifier-extracted, guarded) -> ref_samplesJ.txt
  ref_libm.js                libm attribution harness
  validate.py                validation: mirror, libm, reduction (sympy + intervals),
                             linear identity, Dsep + gradient        -> validation.json
  certify_analytic.py        the analytic closure of section 3       -> analytic_counts.json
  certify_coincident.py      sinc-based reduced forms, certified counts (box path)
  run_coincident.py          driver: box covering of the beta axis (corroboration)
  spot_check.py              outer-interval spot check (corroboration) -> spot_check.json
  certify_separated.py       the 2-D mass-free sweep, torque-row clearing
  certify_separated_wc.py    the same sweep, radial-row clearing (Dred)
  conditioning_measure.py    the sweep's measured conditioning       -> conditioning.json
  run_separated.py           driver for the 2-D sweep over beta windows
  census_cert.py             per-teacher certified returned-family lists:
                             complete sign chart + Krawczyk existence/type;
                             exact separated lists only when paired with the
                             branch-and-bound exclusion cover
  face_bb.py                 the branch-and-bound closure attempt    -> face_bb_attempt.json
  face_bb_signlaw.py         the SIGN-LAW-pruned enumeration (order-only box
                             exclusion licensed by SeparatedCount/*.lean), and
                             the parametrized-teacher mode that certifies the
                             family count constant over a (beta, y) BOX
                                                                    -> signlaw_bb_*.json
  sep_jacobian_audit.py      REGRESSION GUARD for the Krawczyk Jacobian (9.5);
                             run it after any edit to sep_res
  jacobian_fix_delta.py      item-by-item re-run of every stored point under the
                             fixed Jacobian                         -> jacobian_fix_delta.json
  faces_rebuild_fixed.py     the post-fix witness table, provenance kept
                                                                    -> faces_fixed.json
  band_rescan_fixed.py       the post-fix band re-certification      -> bandfix2_*.json
  faces.py                   proposed-label witnesses via locate-then-certify
                                                                    -> faces_0.json
                             (SUPERSEDED by faces_fixed.json, see 9.5)
  wall_pos.py                pointwise mixed-collar witness scan -> wall_pos.json
  locator_probe.json         seed-density doubling probe at the face representatives
  ywall.py / ywall_pos.py    the candidate-transition scans          -> ywallp_0.json
  band_scan.py, map_scan.py, compare_band.py   finite sampled scans; certified
                                               rows remain pointwise
  harness.js, bandJ.js       DOM shim + probe for the SHIPPED widgets.js (read only)
  assemble.py / finalize.py / fill_md.py   builders of certificate.json and of the
                             first-generation CERTIFICATE.md; their per-run inputs
                             (coincident_*.json, separated_*.json, band_*.json,
                             faces_*.json shards) were not retained, so they do not
                             re-run from the committed tree; certificate.json is the
                             committed record of that run
  certificate.json           machine-readable pointwise/box data from the original
                             run plus explicitly proposed arrangement metadata
  band_compare.json          the 158-point band comparison vs widgets.js @ f430fab
  band_comparefix.json       the same 158 points vs widgets.js @ 2201dd0
  CERTIFICATE.md             this file (hand-maintained since the analytic closure)
```
