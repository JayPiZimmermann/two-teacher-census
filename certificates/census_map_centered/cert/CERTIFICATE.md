# Certificate: the boundary set of the centered two-student phase diagram

**Object.** A certified-complete enumeration of the zero sets, in the
teacher-gap / mass-ratio map coordinates `(x = beta in [0,pi], y in (-1,1])`,
of the three mass-free determinants of the centered kernel:

```
  phiCos(t)   = cos t * arcsin(cos t) + |sin t|          (pi-periodic)
  h(t)        = sin t * arcsin(cos t)                    (couplingH)
  slopeAtom(t)= phiCos(t) - 2|sin t| = h'(t)

  Wtau(beta,t) = h(t) sA(t-beta) - h(t-beta) sA(t)          torque Wronskian
  Wpot(beta,t) = phiCos(t) h(t-beta) - phiCos(t-beta) h(t)  potential Wronskian
  Wwgt(beta,t) = |sin t| h(t-beta) - |sin(t-beta)| h(t)     weight Wronskian
```

with `psi = (y+1)pi/2`, `(s0,s1) = (sin psi, cos psi)`, and a root `t` mapped
to `y` through the kernel-vector mass direction `(s0,s1) ~ (-h(t-beta), h(t))`,
folded by the `(-,-)` isometry exactly as `ratioCoord` in `website/widgets.js`.

`phiCos` and `couplingH` are the Lean definitions of
`Planar/SignedAngle/PlanarKernel.lean` -- `phiCos x = orthantPhi (cos x)`
with `orthantPhi rho = rho arcsin(rho) + sqrt(1 - rho^2)`, and
`couplingH D = sin D * arcsin(cos D)` -- and the closed branch forms used by
`widgets.js` and mirrored here, `phiC(t) = (pi/2 - x) cos x + sin x` and
`HC(t) = (pi/2 - x) sin x` with `x = t mod pi`, agree with them for every
`t` (on `[0, pi]`, `arcsin(cos t) = pi/2 - t` and `sqrt(1 - cos^2 t) = sin t`,
and both sides are `pi`-periodic).  The three Wronskians are the Lean
`torqueAtomWronskian`, `potentialAtomWronskian`, `weightAtomWronskian` of
`Planar/GeneralTeachers/PhaseBoundary.lean` / `PotentialBoundary.lean`,
definition for definition.

**Scope caveat.** This certificate is about the zero sets of these three
determinants only. It does **not** claim that they exhaust the census
boundary of the classification; that containment is a separate (Lean) matter.

**Precision.** All reported enclosures come from `mpmath.iv` (interval
arithmetic with directed rounding) at **200 bits** of working precision
(mpmath 1.3.0). `float64` is used only inside the widgets.js mirror of section 1
and in pre-scans that guess where roots are; no reported enclosure passes
through it.

## 1. Validation against `website/widgets.js`

The reference values are produced by `ref_widgets.js`: `node` evaluating the
**literal function sources** of `website/widgets.js` (`PI`, `mod`, `phiC`,
`HC`, `dHC`, `slopeAtom`), extracted by identifier with a guard that asserts
each extracted function still contains the expected expression, and `eval`-ed
verbatim.  Last re-run against `widgets.js` at commit `2201dd0` (whose
centered kernel functions are byte-identical to those first validated at
`f430fab`).
Sample count: **25000** quasi-random `(t, beta)` pairs, `t` in `[-2pi, 3pi]`,
`beta` in `[0, pi]`.

| quantity | max abs deviation, Python mirror vs widgets.js |
|---|---|
| `phiC` | 4.441e-16 |
| `HC` | 1.110e-16 |
| `dHC` | 4.441e-16 |
| `slopeAtom` | 4.441e-16 |
| `Wtau` | 3.331e-16 |
| `Wpot` | 4.441e-16 |
| `Wwgt` | 2.220e-16 |

These deviations are **1-2 ulp and are entirely attributable to libm**
(`ref_libm.js` + section (a2) of `validate.py`): V8's `Math.sin`/`Math.cos`
and glibc's differ by up to 1.1102e-16 (measured over 6000 arguments), and feeding
node's *own* `sin`/`cos` values into the Python formulas reproduces every
one of the seven quantities **bit-identically (max deviation 0.0)**, so
the reimplementation is formula-identical to the widget.

`h'(t) = slopeAtom(t)` by central differences (central difference, mpmath prec=300, step 1e-25, kink lattice avoided by 1e-3):
max absolute error **4.9998537e-51**, attained at `t = 11.0015`.

## 2. An exact reduction that removes the kinks

`|sin|` and the pi-periodic folding put kinks at `t in pi*Z` and
`t - beta in pi*Z`. An interval derivative test is invalid across a kink, so
the `t`-period `[0,pi]` is **split at the kink lattice** first. For
`beta in (0,pi)` the split points are exactly

```
    t = 0 ,   t = beta ,   t = pi
```

giving two smooth pieces

```
    piece B : t in [0, beta]   fold(t) = t, fold(t-beta) = t-beta+pi
    piece A : t in [beta, pi]  fold(t) = t, fold(t-beta) = t-beta
```

On each closed piece the widgets.js *closed branch forms* are the correct
values (the one-sided form is still correct at the branch endpoint), so the
piece is analytic and the interval derivative test is legitimate. In the code
this is `centered.dets_branch(B, T, kt, ku)`, which takes the branch indices
`(kt, ku) = (0,-1)` on piece B and `(0,0)` on piece A rather than inferring
them; the generic multi-branch evaluator `centered.atoms` (hull over all
branches met, with a runtime assertion that the branch pieces tile the input)
is used only where a value, not a derivative, is needed.

**Worked example.** At `beta = pi/2`, the split points are `t = 0, pi/2, pi`.
On piece B the enclosure of `Wtau` over `t in [0, pi/2]` contains 0 *and* the
enclosure of `d/dt Wtau` contains 0, so the naive test cannot conclude. Split
at the kink and work on the closed piece: the reduced form below shows
`Phi(pi/2, w) > 0` for `|w| < pi/4` and `Phi(pi/2, +-pi/4) = 0` with
`Phi_w(pi/2, +-pi/4) = 0` as well -- the root sits *exactly on* the kink and
is a double root there. Had the domain not been split, the derivative
enclosure across `t = 0` would have been taken over a set on which `Wtau` is
not differentiable, and the resulting 'monotone, hence at most one root'
conclusion would have been unsound.

Substituting the branch forms and putting, on each piece,

```
   piece A:  b = beta   , sigma = +1 , w = t - (beta+pi)/2
   piece B:  b = pi-beta, sigma = -1 , w = t - beta/2
   W = (pi - b)/2 ,  |w| <= W

   Q (b,w) = (w^2 - b^2/4) sin b
   Om(b,w) = (b/2)(cos 2w + cos b) = b cos(w + b/2) cos(w - b/2)

   Wtau = sigma * ( Q + Om) =: sigma * Phi(b,w)
   Wpot = sigma * (-Q + Om) =: sigma * Psi(b,w)
   Wwgt = sigma *       Om
```

`Phi`, `Psi`, `Om` are **analytic in `(b,w)` on the whole piece and even in
`w`**: the reduction has absorbed every kink into the piece boundary
`|w| = W`. This is verified two ways in `validate.py`: symbolically with
`sympy` on both pieces ([True, True, True] / [True, True, True], all six identities), and on **5782 interval
samples** where the enclosure of (folded mirror - reduced form) contains 0
with radius at most 7.84e-59.

### 2.1 An exact linear relation

Because `Phi + Psi = 2 Om` identically,

```
    Wtau(beta,t) + Wpot(beta,t) = 2 Wwgt(beta,t)      for all beta, t
```

proved symbolically (`sympy`, exact) and verified in interval arithmetic with
residual radius at most 1.02e-58. **The three determinants span only a
2-dimensional space**; any two of them determine the third.

## 3. Structure lemmas (the analytic backbone)

Endpoint closed forms, verified at 399 grid values of `b` in interval
arithmetic (max residual 3.48e-59):

```
   Phi(b,0)   = b cos(b/2) [cos(b/2) - (b/2) sin(b/2)]
   Phi(b,W)   = (pi/4)(pi - 2b) sin b
   Psi(b,0)   = b cos(b/2) [cos(b/2) + (b/2) sin(b/2)]   > 0 on (0,pi)
   Psi(b,W)   = (pi/2)(b - pi/2) sin b
   Om(b,W)    = 0
```

**L5.** Psi_w = -(2w sin b + b sin 2w) < 0 for 0 < w <= W and 0 < b < pi, so Psi is strictly decreasing in w on [0,W].

**L6.** Om = b cos(w+b/2) cos(w-b/2); for |w| <= W both arguments lie in [b-pi/2, pi/2] subset [-pi/2,pi/2], hence Om >= 0 with equality exactly at |w| = W.

**L7.** for b > pi/2, D := -Phi_w = b sin 2w - 2w sin b satisfies D(0)=0, D_ww = -4b sin 2w <= 0 (concave on 0<=2w<=pi-b) and D(W) = (2b-pi) sin b > 0, hence D(w) >= (w/W) D(W) > 0: Phi is strictly decreasing in w on [0,W].

**L8.** for 0 < b < pi/2, Phi > 0 on [0,W].  For w in [b/2,W], Phi = Q+Om with Q,Om >= 0 and never both zero.  For w in [0,b/2] put p = b/2+w, q = b/2-w >= 0, p+q = b <= pi/2; then Phi = (p+q)cos p cos q - pq sin(p+q) and, dividing by cos p cos q > 0, the claim is p+q > pq(tan p + tan q). Since q <= pi/2-p one has pq tan p <= p(pi/2-p)tan p < p, the last step being (pi/2-p)tan p < 1, i.e. s < tan s at s = pi/2-p; symmetrically pq tan q < q.  Adding gives the claim strictly unless p=q=0.

Interval grid checks of these hypotheses: `{'L5': 7960, 'L6': 8159, 'L7': 3960, 'L8': 4059}`.

## 4. Certified root counts as a function of beta

Roots are counted in `t` modulo `pi` (all three determinants are
`pi`-periodic in `t`). The counting is done on `w in [0,W]` and doubled by
the exact evenness in `w`.

| beta | Wtau | Wpot | Wwgt |
|---|---|---|---|
| `beta = 0` | identically zero | identically zero | identically zero |
| `(0, beta*_1)` | 0 | 2 | 2 |
| `beta = beta*_1` | **1 (double root)** | 2 | 2 |
| `(beta*_1, pi/2)` | 2 | 2 | 2 |
| `beta = pi/2` | 2 (both on the kink lattice) | 2 (same two points) | 2 (same two points) |
| `(pi/2, beta*_2)` | 2 | 2 | 2 |
| `beta = beta*_2` | **1 (double root)** | 2 | 2 |
| `(beta*_2, pi)` | 0 | 2 | 2 |
| `beta = pi` | identically zero | identically zero | identically zero |

Only **one** count changes anywhere: that of `Wtau`, and only at
`beta = beta*_1` and `beta = beta*_2`, in both cases through a **double root**
(a common zero of `Phi(b,.)` and `Phi_w(b,.)` at `w = 0`, i.e. `Phi_w(b,0) = 0`
holds identically and `Phi(b*,0) = 0`). `Wpot` and `Wwgt` have constant count
2 on all of `(0,pi)`.

Certified coverage of the beta axis:

| beta range | boxes | measure covered | uncovered |
|---|---|---|---|
| `(0, beta*_1)` | 80 | 1.420925475550880 | 5.396e-14 |
| `(beta*_1, pi/2)` | 94 | 0.149870851243681 | 1.817e-13 |
| `(pi/2, beta*_2)` | 94 | 0.149870851243681 | 1.817e-13 |
| `(beta*_2, pi)` | 80 | 1.420925475550880 | 5.396e-14 |

The uncovered part consists only of collars around the critical betas
`0, beta*_1, pi/2, beta*_2, pi`, each of width below `1e-13`; those five
betas themselves are handled exactly (section 5). See section 9.

## 5. The critical betas

### `beta = 0`
- **nature**: degenerate column
- **statement**: t - beta = t, so Wtau = Wpot = Wwgt = 0 identically in t; the entire column beta = 0 belongs to all three zero sets.
- **count**: uncountable (identically zero)

### `beta = beta*_1`
- **nature**: double root of the torque Wronskian (lens endpoint)
- **beta_enclosure**: ['1.420925475551033713494856535', '1.420925475551033713494856535']
- **double_root_t**: t = beta/2   (piece B, w = 0)
- **why**: Phi(pi-beta,0) = 0 and Phi_w(b,0) = 0 identically, so w = 0 is a double root; the Wtau count jumps 0 -> 2 here.
- **Phi_and_Phi_w_enclosures**: [['-1.7424443e-59', '1.617984e-59'], ['0.0', '0.0']]
- **y_of_the_double_root**: -1/2 exactly (h(beta/2) and -h(-beta/2) are equal)

### `beta = pi/2`
- **nature**: total degeneracy / orthogonal teacher
- **statement**: b = pi/2 on BOTH pieces and W = pi/4.  Phi(pi/2,pi/4) = Psi(pi/2,pi/4) = Om(pi/2,pi/4) = 0 and Phi_w(pi/2,pi/4) = 0, so all three determinants vanish exactly on the same set {t = 0, t = pi/2} (mod pi) -- both points lying ON the kink lattice -- and there is no interior root (L5/L7 with Phi(.,0), Psi(.,0) > 0).
- **values_at_(pi/2,pi/4)**: [['-5.9543425e-60', '3.556407e-60'], ['-5.332041e-60', '4.1787086e-60'], ['-2.8428349e-60', '1.0672009e-60']]
- **Phi_w_at_(pi/2,pi/4)**: ['-3.7338092e-60', '3.7338092e-60']
- **count**: {'Wtau': 2, 'Wpot': 2, 'Wwgt': 2}
- **mass_direction**: at both roots h(t) = h(t-beta) = 0, so the mass direction is 0/0; the curves extend by the limits y = c-1, -c (torque) and y = 1-c, c (potential).

### `beta = beta*_2`
- **nature**: double root of the torque Wronskian (lens endpoint)
- **beta_enclosure**: ['1.72066717803875952496778684828', '1.72066717803875952496778684828']
- **double_root_t**: t = (beta+pi)/2   (piece A, w = 0)
- **why**: Phi(beta,0) = 0 and Phi_w(b,0) = 0 identically; the Wtau count jumps 2 -> 0 here.
- **Phi_and_Phi_w_enclosures**: [['-1.0579126e-59', '9.9568244e-60'], ['0.0', '0.0']]
- **y_of_the_double_root**: -1/2 exactly

### `beta = pi`
- **nature**: degenerate column (identified with beta = 0 by pi-periodicity)
- **statement**: t - beta = t (mod pi); all three determinants vanish identically.

## 6. beta*: independent verification of the lens claim

The lens endpoints are the double-root locus `Phi(b,0) = 0`. Since
`Phi(b,0) = b cos(b/2) [cos(b/2) - (b/2) sin(b/2)]` and `cos(b/2) > 0` on
`(0,pi)`, this is `cos(b/2) = (b/2) sin(b/2), i.e. (b/2) tan(b/2) = 1`.

Uniqueness in the bracket `['1.70', '1.75']`: f'(b) = -sin(b/2) - (b/4) cos(b/2) certified strictly negative on [1.70,1.75], and f(1.70) > 0 > f(1.75).

```
   b*      in 1.72066717803875952496778684828
   beta*_1 = pi - b*  in 1.420925475551033713494856535
   beta*_2 = b*       in 1.72066717803875952496778684828
   beta*_1 + beta*_2 - pi  in [-1.2446e-59, 1.2446e-59]   (so the pair is symmetric about pi/2)
```

**The previous lane's characterisation survives.** With
`u* = beta*_1 / 2 in 0.710462737775516856747428267502`, the residual `tan u* - (pi/2 - u*)` is enclosed in
`[-8.71222e-60, 8.71222e-60]`, i.e. certified to contain 0.
The two characterisations are *equivalent*, not merely numerically close:
`u* = pi/2 - b*/2`, so `tan u* = cot(b*/2) = cos(b*/2)/sin(b*/2)`, and
`tan u* = pi/2 - u* = b*/2` is exactly `cos(b*/2) = (b*/2) sin(b*/2)`.

The reported value `beta* = 1.420925475551` agrees with the certified
enclosure to all quoted digits.

**Lens claim: SURVIVES, as a single closed component.** For every
`beta in (beta*_1, beta*_2)` the certified `Wtau` count is exactly 2; at both
endpoints it is exactly 1 (a double root); outside `[beta*_1, beta*_2]` it is
0. The two roots are `w = +-w_r(beta)` with `w_r` the unique zero of a
function certified strictly monotone in `w` near it, so `w_r` is continuous
(implicit function theorem) and the zero set is one closed curve. The two
branches are exact mirror images: `t -> beta - t (mod pi)` maps roots to roots
and sends the mass ratio `rho` to `1/rho`, hence `y` to `-1-y`, so the lens is
**symmetric in `y` about `-1/2`** and both endpoints lie exactly on `y = -1/2`.
This identity was checked on every certified curve box: `y+ + y- = -1`: **True**.

`y`-hull over the certified curve boxes: `[-0.63958222380843369103, -0.36041777619655092697]`
(beta measure covered 0.299691574175 of 0.299741702488).
The extremes `y = c-1` and `y = -c` with `c = (2/pi) arctan(2/pi) in 0.360907073228108373304919068364`
are the `beta -> pi/2` limits and are not attained.

## 7. The 1-dimensional pieces

**1D-1. torque Wronskian lens, upper branch**
  - equation: Wtau(beta,t) = 0, the root whose mass ratio has rho > 1 (y > -1/2)
  - beta_range: [beta*_1, beta*_2] = [1.420925475551033713495, 1.720667178038759524968]
  - parameterisation: the two roots are w = +-w_r(beta) with w_r the unique zero of Phi(b,.) in (0,(pi-b)/2); on piece B (beta < pi/2, b = pi-beta, t = beta/2 + w) this branch is w = +w_r, on piece A (beta > pi/2, b = beta, t = (beta+pi)/2 + w) it is w = -w_r
  - y_range: [ -1/2 , -c )  with c = 0.3609070732281083733
  - certified_boxes: certificate.json: curves.torque_lens.certified_boxes_sample (full list: certificate_full.json.gz)
  - note: y = -1/2 exactly at the two endpoints beta*_1, beta*_2; y -> -c as beta -> pi/2 (limit, not attained)

**1D-2. torque Wronskian lens, lower branch**
  - equation: Wtau(beta,t) = 0, the root with rho < 1 (y < -1/2)
  - beta_range: same as 1D-1
  - parameterisation: exact mirror of 1D-1 under t -> beta - t (mod pi), which sends the mass ratio rho to 1/rho and hence y to -1-y
  - y_range: ( c-1 , -1/2 ]
  - certified_boxes: certificate.json: curves.torque_lens.certified_boxes_sample (full list: certificate_full.json.gz)

**1D-3. potential Wronskian curve, lower branch**
  - equation: Wpot(beta,t) = 0, the root with y < 1/2
  - beta_range: (0, pi)  (all of it)
  - parameterisation: the two roots are w = +-w_r(beta) with w_r the unique zero of Psi(b,.) in (0,(pi-b)/2); piece A for beta < pi/2 (b = beta, t = (beta+pi)/2 + w, this branch is w = -w_r), piece B for beta > pi/2 (b = pi-beta, t = beta/2 + w, this branch is w = +w_r)
  - y_range: ( c , 1/2 )  -- entirely inside the MIXED-sign sector y > 0
  - certified_boxes: certificate.json: curves.potential_curve.certified_boxes_sample (full list: certificate_full.json.gz)

**1D-4. potential Wronskian curve, upper branch**
  - equation: Wpot(beta,t) = 0, the root with y > 1/2
  - beta_range: (0, pi)
  - parameterisation: mirror of 1D-3 under t -> beta - t (mod pi); y -> 1-y
  - y_range: ( 1/2 , 1-c )
  - certified_boxes: certificate.json: curves.potential_curve.certified_boxes_sample (full list: certificate_full.json.gz)

**1D-5. weight Wronskian image, root t = 0 (mod pi)**
  - equation: Wwgt(beta,t) = 0 at t = 0; there h(t) = 0, so s1 = 0
  - beta_range: (0, pi)
  - image_in_map: the horizontal line y = 0
  - note: coincides with the degenerate stratum 1D-9 (second teacher massless).  Here s1 = h(0) = 0 exactly and s0 != 0, so atan2 gives +-pi/2 and the fold gives y = 0; the per-beta table reports an enclosure of width < 1e-27 around 0.

**1D-6. weight Wronskian image, root t = beta (mod pi)**
  - equation: Wwgt(beta,t) = 0 at t = beta; there h(t-beta) = 0, so s0 = 0
  - beta_range: (0, pi)
  - image_in_map: the horizontal line y = +1 (identified with y = -1)
  - note: coincides with the degenerate stratum 1D-10.  In the per-beta table the y field of this root is reported as null: s0 = -h(0) = 0 EXACTLY, so atan2(s0,s1) sits on the seam of the (-,-) fold and the interval routine correctly declines to pick a branch.  The value y = +1 follows from the exact vanishing of s0 together with ratioCoord's convention a <= 0 -> a + pi, not from a numerical evaluation.

**1D-7. degenerate column beta = 0**
  - equation: beta = 0
  - content: all three determinants vanish identically in t; coincident teacher directions

**1D-8. degenerate column beta = pi/2 (orthogonal teacher)**
  - equation: beta = pi/2
  - content: all three determinants have the SAME zero set {t = 0, t = pi/2} (mod pi), both on the kink lattice; the mass direction degenerates to 0/0 at both

**1D-9. degenerate stratum y = 0**
  - equation: s1 = 0 -- the second teacher atom is massless

**1D-10. degenerate stratum y = +1 (= y = -1)**
  - equation: s0 = 0 -- the first teacher atom is massless; psi = 0 and psi = pi are the same projective line, so the top and bottom edges of the box are one stratum

**1D-11. degenerate column beta = pi**
  - equation: beta = pi
  - content: identified with beta = 0 by the pi-periodicity of the centered kernel; all three determinants vanish identically

### Non-intersections (certified)

- The torque lens (1D-1/1D-2) and the potential curve (1D-3/1D-4) never meet in the map: a common zero of Wtau and Wpot forces Q = Om = 0, and Om = 0 only at |w| = W, where Q = (pi/4)(pi-2b) sin b, so b = pi/2, i.e. beta = pi/2 and t in {0, pi/2}.  At that (beta,t) the mass direction is 0/0 and the two curves approach DIFFERENT limits (y = -c, c-1 versus y = c, 1-c), so their images are disjoint.
- The torque lens stays inside y in (c-1, -c) subset (-1,0), hence never meets y = 0 or y = +-1.
- The potential curve stays inside y in (c, 1-c) subset (0,1), hence never meets y = 0 or y = +-1.

## 8. The 0-dimensional pieces

| id | name | beta | y | what makes it special |
|---|---|---|---|---|
| 0D-1 | left lens endpoint / torque double root | `1.420925475551033713495` | `-0.5` | double root of Wtau at t = beta/2 (w = 0): Phi(pi-beta,0) = 0 and Phi_w(.,0) = 0.  The Wtau root count jumps 0 -> 2.  Meeting point of 1D-1 and 1D-2.  beta*_1 = 2u* with tan u* = pi/2 - u*. |
| 0D-2 | right lens endpoint / torque double root | `1.720667178038759524968` | `-0.5` | double root of Wtau at t = (beta+pi)/2 (w = 0).  Count jumps 2 -> 0.  beta*_2 = pi - beta*_1 = b*, (b*/2) tan(b*/2) = 1. |
| 0D-3 | torque lens meets the column beta = pi/2 (upper) | `1.570796326794896619231` | `-0.3609070732281083733` | the Wtau root reaches the kink t = beta (equivalently t = 0); the mass direction is 0/0 there and the curve extends by the limit y = -c, c = (2/pi) arctan(2/pi).  y-extremum of 1D-1. |
| 0D-4 | torque lens meets the column beta = pi/2 (lower) | `1.570796326794896619231` | `-0.6390929267718916267` | mirror of 0D-3 under y -> -1-y; y-extremum of 1D-2. |
| 0D-5 | potential curve meets the column beta = pi/2 (lower) | `1.570796326794896619231` | `0.3609070732281083733` | the Wpot root reaches the kink lattice; limit y = c.  y-extremum of 1D-3. |
| 0D-6 | potential curve meets the column beta = pi/2 (upper) | `1.570796326794896619231` | `0.6390929267718916267` | mirror under y -> 1-y; y-extremum of 1D-4. |
| 0D-7 | potential curve endpoint on the column beta = 0 | `0` | `0.5` | as beta -> 0 the two Wpot roots tend to t = pi/2 +- w_D with w_D the Dottie number (cos w = w), and the mass ratio tends to rho = -1, i.e. y = 1/2.  1D-3 and 1D-4 meet here. |
| 0D-8 | potential curve endpoint on the column beta = pi | `3.141592653589793238463` | `0.5` | mirror of 0D-7 under beta -> pi - beta. |
| 0D-9 | beta = 0 meets y = 0 | `0` | `0` | corner of the degenerate strata 1D-7 and 1D-9 (also the endpoint of the weight-Wronskian image 1D-5) |
| 0D-10 | beta = 0 meets y = +1 (= -1) | `0` | `1` | corner of 1D-7 and 1D-10 (endpoint of 1D-6) |
| 0D-11 | beta = pi/2 meets y = 0 | `1.570796326794896619231` | `0` | crossing of the degenerate column 1D-8 with 1D-9/1D-5 |
| 0D-12 | beta = pi/2 meets y = +1 (= -1) | `1.570796326794896619231` | `1` | crossing of 1D-8 with 1D-10/1D-6 |
| 0D-13 | beta = pi meets y = 0 | `3.141592653589793238463` | `0` | corner of 1D-11 and 1D-9 |
| 0D-14 | beta = pi meets y = +1 (= -1) | `3.141592653589793238463` | `1` | corner of 1D-11 and 1D-10 |

Auxiliary constants appearing above:

```
   u*   in 0.710462737775516856747428267502      (tan u* = pi/2 - u*, beta*_1 = 2u*)
   b*   in 1.72066717803875952496778684828      ((b*/2) tan(b*/2) = 1, beta*_2 = b*)
   w_D  in 0.739085133215160641655312087674      (Dottie number, cos w = w)
   c    in 0.360907073228108373304919068364      (c = (2/pi) arctan(2/pi))
```

## 9. Unresolved set

Total entries: **55440**.

### 9.1 beta collars where a root COUNT could not be certified

| beta range | box | width |
|---|---|---|
| `(0, beta*_1)` | `[1.420925475550979866654, 1.420925475551033713495]` | 5.385e-14 |
| `(beta*_1, pi/2)` | `[1.420925475551033713495, 1.420925475551124584658]` | 9.087e-14 |
| `(beta*_1, pi/2)` | `[1.570796326794805748068, 1.570796326794896619231]` | 9.087e-14 |
| `(pi/2, beta*_2)` | `[1.570796326794896619231, 1.570796326794987490395]` | 9.087e-14 |
| `(pi/2, beta*_2)` | `[1.720667178038668653804, 1.720667178038759524968]` | 9.087e-14 |
| `(beta*_2, pi)` | `[1.720667178038759524968, 1.720667178038813371809]` | 5.385e-14 |

Total uncovered beta measure: **4.712e-13**.
Each collar sits around one of the critical betas `beta*_1, pi/2, beta*_2`,
which are themselves treated exactly in section 5. This is the ONLY part of
the beta axis on which a root count is not certified.

### 9.2 curve boxes whose y enclosure is not tight

These are boxes on which the ROOT is certified to exist and be unique, but
the resulting enclosure of `y` is wider than the reporting tolerance `1e-3`,
or the mass direction `(-h(t-beta), h(t))` is not sign-determined on the box.
They are certified, just not tight; they are excluded from the `y`-hulls of
section 6.

| curve | wide boxes | beta measure | max box width | reasons |
|---|---|---|---|---|
| torque lens | 22864 | 0.000050 | 8.980e-09 | count: 4, mass_direction_undetermined: 40, y_enclosure_wider_than_0.001: 22820 |
| potential curve | 32570 | 0.000636 | 9.412e-08 | count: 2, mass_direction_undetermined: 58, y_enclosure_wider_than_0.001: 32510 |

The `mass_direction_undetermined` boxes (98 in total) are exactly the ones
adjacent to `beta = pi/2`, where the defining root reaches the kink lattice
and `h(t) = h(t-beta) = 0`; the `count` boxes (6) are the collars of 9.1.

### 9.3 the generic root test

The generic exclusion/monotonicity/bisection test of section 5 left **no**
undecided box that is not adjacent to a kink point: every box it could not
decide abuts `t = 0`, `t = beta` or `t = pi`, where the determinant is
certified to vanish (Wwgt) or certified nonzero (Wtau, Wpot) by the closed
forms of section 3. Those boxes are reported per beta in
`per_beta_table[*][det].kink_adjacent_boxes_explained`.

## 10. Honesty notes

- Every count in section 4 is certified **for all beta in the stated box**, not
  sampled. Where a count could not be certified, the box is listed in section 9
  rather than rounded to the nearest plausible integer.
- The `y`-extremes `c-1, -c` (torque) and `c, 1-c` (potential) are **limits**
  at `beta = pi/2`, derived analytically; at `beta = pi/2` itself the defining
  root sits on the kink lattice with `h(t) = h(t-beta) = 0`, so the mass
  direction is a genuine `0/0` and the point is degenerate.
- The claim that these three determinants exhaust the census boundary is **not**
  made here.
- Validity domain: `beta in [0, pi]`, `t` in one `pi`-period, centered kernel
  only. Precision 200 bits throughout.

## Files

```
  centered.py                kernel module: float mirror of widgets.js + mpmath.iv interval forms
  ref_widgets.js             reference generator: literal widgets.js functions (identifier-extracted, guarded) -> ref_samples.txt
  ref_libm.js                libm attribution harness (node's own sin/cos per sample)
  validate.py                validation: widgets.js mirror, libm attribution, reduction proof (sympy + intervals), linear identity
  validation.json            the latest validation run (regenerate: node ref_widgets.js && python3 validate.py)
  certify.py                 the certification
  compact.py                 downsamples the raw output into certificate.json
  make_md.py                 regenerates CERTIFICATE.md from certificate.json + validation.json
  certificate.json           the compact machine-readable certificate (this file)
  certificate_full.json.gz   the complete run, every certified and wide box (not in version control; certify.py regenerates it)
  CERTIFICATE.md             human-readable appendix
  run_d16_shipped.log        console log of the certification run
```

## Provenance

```
  MAX_DEPTH_Y                     16
  console log                     run_d16_shipped.log
  sha256 certificate.json         7737a899aae6b4fa9263c5206395743d484649204903063b44d82425ec55422f
  sha256 certificate_full.json.gz 230e7f4a52ff3882bce17b04aa3c9271ce5e0be4568c8860f2d5b40e04040941
```

The depth of the shipped run is recorded in `cert/PROVENANCE.json` rather
than in `certificate.json` itself: this certificate predates the change
that stamps the run parameters into every emitted certificate (`certify.py`
now writes them under `run_parameters`), and nothing was injected into the
shipped JSON after the fact.
