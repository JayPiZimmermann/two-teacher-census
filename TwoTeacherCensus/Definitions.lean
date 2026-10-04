import Mathlib

/-!
# The planar two-atom setting: definitions

Vocabulary of the two-student / two-teacher census, in the notation of the
paper repository (`Removing-spurious-minima-for-planar-features-by-skip-connections/
Definitions.lean`).  Everything here is self-contained over Mathlib; the kernel
definitions `orthantPhi`, `reluPhi`, `phiCos`, `phiCosJ`, `couplingH`,
`couplingHJ`, `kernelMassLoss`, `kernelTeacherSelfEnergy`, `kernelExcessLoss`
are copied verbatim from that file so that the two developments share names.

Argument order for the two-atom objects follows the manuscript tables and the
long tree: teacher directions `β`, teacher masses `s`, student masses `c`,
student directions `θ` (all `Fin 2 → ℝ`).  In the general-width kernel
functionals copied from the paper repository the order is the paper's
(`student masses, teacher angles, teacher masses, student angles`).

Conventions.
* `Φ q` is the angular kernel `Φ_q`, `Hq q = -Φ_q'` the coupling, `κ q = Φ_q 0`.
* The population loss `L q β s c θ` is the finite kernel form
  `(1/(2π)) · kernelExcessLoss (Φ q) c β s θ`; the Gaussian-integral identity
  that identifies it with the literal planar population loss is proved in the
  paper repository (`eq_loss_quadratic`) and is *not* re-derived here.
* Criticality is `fderiv = 0` on the product `(Fin 2 → ℝ) × (Fin 2 → ℝ)`;
  local minima are Mathlib's `IsLocalMin` in all four coordinates, with signed
  masses allowed in every neighbourhood.
* A **topological saddle** is a critical point that is neither a local minimum
  nor a local maximum.
-/

noncomputable section

open Real Set Filter Topology Classical
open scoped BigOperators

namespace TwoTeacherCensus

/-! ## The two features -/

/-- The two hidden-unit features: centered `C(z) = |z|/2` and plain ReLU. -/
inductive Model where
  | centered
  | plainRelu
  deriving DecidableEq

/-! ## Scalar kernels (verbatim from the paper repository) -/

def orthantPhi (ρ : ℝ) : ℝ := ρ * arcsin ρ + sqrt (1 - ρ ^ 2)
def reluPhi (ρ : ℝ) : ℝ := orthantPhi ρ + (π / 2) * ρ
def phiCos (x : ℝ) : ℝ := orthantPhi (cos x)
def phiCosJ (x : ℝ) : ℝ := reluPhi (cos x)
def couplingH (x : ℝ) : ℝ := sin x * arcsin (cos x)
def couplingHJ (x : ℝ) : ℝ := couplingH x + (π / 2) * sin x

/-- The angular kernel `Φ_q` of the manuscript. -/
def Φ : Model → ℝ → ℝ
  | .centered => phiCos
  | .plainRelu => phiCosJ

/-- The angular coupling `H_q = -Φ_q'`. -/
def Hq : Model → ℝ → ℝ
  | .centered => couplingH
  | .plainRelu => couplingHJ

/-- The diagonal kernel value `κ_q = Φ_q(0)`: `π/2` and `π`. -/
def κ : Model → ℝ
  | .centered => π / 2
  | .plainRelu => π

/-- The kernel period: `π` (projective directions) and `2π` (oriented rays). -/
def period : Model → ℝ
  | .centered => π
  | .plainRelu => 2 * π

/-- The slope of the coupling off the kink lattice, `H_q' = Φ_q - 2|sin|`
(the kernel equation `Φ_q'' + Φ_q = 2|sin|`). -/
def Hslope (q : Model) (x : ℝ) : ℝ := Φ q x - 2 * |sin x|

/-! ## Finite kernel energies (verbatim from the paper repository) -/

def kernelMassLoss {n m : ℕ} (kf : ℝ → ℝ) (s : Fin n → ℝ) (β : Fin m → ℝ)
    (t : Fin m → ℝ) (θ : Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ i, ∑ j, s i * s j * kf (θ i - θ j) -
    ∑ i, ∑ k, s i * t k * kf (θ i - β k)

def kernelTeacherSelfEnergy {m : ℕ} (kf : ℝ → ℝ) (β : Fin m → ℝ)
    (t : Fin m → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∑ k, ∑ l, t k * t l * kf (β k - β l)

def kernelExcessLoss {n m : ℕ} (kf : ℝ → ℝ) (s : Fin n → ℝ) (β : Fin m → ℝ)
    (t : Fin m → ℝ) (θ : Fin n → ℝ) : ℝ :=
  kernelMassLoss kf s β t θ + kernelTeacherSelfEnergy kf β t

/-! ## The two population losses in planar kernel form -/

/-- The population loss `L_q` of a width-`n` planar student `(c, θ)` against a
width-`m` planar teacher `(s, β)`, in the exact finite kernel form
`L_q = (1/(2π)) [½ Σ c_i c_j Φ_q(θ_i-θ_j) - Σ c_i s_k Φ_q(θ_i-β_k) + ½ Σ s_k s_l Φ_q(β_k-β_l)]`.
It is nonnegative and vanishes exactly at the exact fits. -/
def L {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : ℝ :=
  (1 / (2 * π)) * kernelExcessLoss (Φ q) c β s θ

/-- The centered population loss `L_C`. -/
abbrev L_C {n m : ℕ} (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : ℝ := L .centered β s c θ

/-- The plain-ReLU population loss `L_R`. -/
abbrev L_R {n m : ℕ} (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : ℝ := L .plainRelu β s c θ

/-- The loss as a function of the joint (mass, angle) parameter. -/
def lossPair {n m : ℕ} (q : Model) (β s : Fin m → ℝ)
    (p : (Fin n → ℝ) × (Fin n → ℝ)) : ℝ :=
  L q β s p.1 p.2

/-! ## Variational vocabulary -/

/-- Criticality in all `2n` coordinates: the Fréchet derivative vanishes. -/
def IsCritical {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  fderiv ℝ (lossPair q β s) (c, θ) = 0

/-- Compatibility spellings matching the long tree. -/
abbrev IsCriticalL_C {n m : ℕ} (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsCritical .centered β s c θ
abbrev IsCriticalL_R {n m : ℕ} (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsCritical .plainRelu β s c θ

/-- Local minimality in all coordinates, signed masses allowed. -/
def IsLocalMinL {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsLocalMin (lossPair q β s) (c, θ)

/-- Local maximality in all coordinates. -/
def IsLocalMaxL {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsLocalMax (lossPair q β s) (c, θ)

/-- Global minimality at fixed student width. -/
def IsGlobalMin {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  ∀ c' θ' : Fin n → ℝ, L q β s c θ ≤ L q β s c' θ'

/-- A nonglobal (spurious) local minimum. -/
def IsSpuriousMin {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsLocalMinL q β s c θ ∧ ¬ IsGlobalMin q β s c θ

/-- A topological saddle: critical, not a local minimum, not a local maximum. -/
def IsTopologicalSaddle {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsCritical q β s c θ ∧ ¬ IsLocalMinL q β s c θ ∧ ¬ IsLocalMaxL q β s c θ

/-- The exhaustive three-way label of the manuscript
(`eq-classification-three-labels`). -/
def VariationalType {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : Prop :=
  IsGlobalMin q β s c θ ∨ IsSpuriousMin q β s c θ ∨ IsTopologicalSaddle q β s c θ

/-! ## Two-atom teacher data (`eq-two-student-teacher-data`, `eq-classification-tau`) -/

/-- The teacher potential `P_q(t) = Σ_k s_k Φ_q(t - β_k)`. -/
def P (q : Model) (β s : Fin 2 → ℝ) (t : ℝ) : ℝ := ∑ k, s k * Φ q (t - β k)

/-- The angular derivative (torque) `A_q = -P_q' = Σ_k s_k H_q(t - β_k)`. -/
def A (q : Model) (β s : Fin 2 → ℝ) (t : ℝ) : ℝ := ∑ k, s k * Hq q (t - β k)

/-- The transverse weight `W(t) = Σ_k s_k |sin(t - β_k)|` (feature-independent). -/
def W (β s : Fin 2 → ℝ) (t : ℝ) : ℝ := ∑ k, s k * |sin (t - β k)|

/-- The negative potential curvature `τ_q = P_q - 2W` (`P_q'' = -τ_q` off the kinks). -/
def τ (q : Model) (β s : Fin 2 → ℝ) (t : ℝ) : ℝ := P q β s t - 2 * W β s t

/-- The slope `W'` of the transverse weight on a smooth cell,
`Σ_k s_k sign(sin(t-β_k)) cos(t-β_k)`. -/
def Wslope (β s : Fin 2 → ℝ) (t : ℝ) : ℝ :=
  ∑ k, s k * Real.sign (sin (t - β k)) * cos (t - β k)

/-- `t` is off the kink lattice `{β_k} + πℤ`. -/
def OffKink (β : Fin 2 → ℝ) (t : ℝ) : Prop := ∀ k, sin (t - β k) ≠ 0

/-- The radial determinant `Δ_q(D) = κ_q² - Φ_q(D)²` (`eq-classification-radial-determinant`). -/
def Δ (q : Model) (D : ℝ) : ℝ := κ q ^ 2 - Φ q D ^ 2

/-- The residual field `F_q(γ) = Σ_i c_i Φ_q(γ-θ_i) - P_q(γ)` of a two-student
network; its value at `θ_i` is the radial row `R_i`. -/
def residual (q : Model) (β s c θ : Fin 2 → ℝ) (γ : ℝ) : ℝ :=
  ∑ i, c i * Φ q (γ - θ i) - P q β s γ

/-- Cramér solutions of the two radial rows for a separated pair
(`eq-separated-cramer`). -/
def cramer0 (q : Model) (β s θ : Fin 2 → ℝ) : ℝ :=
  (κ q * P q β s (θ 0) - Φ q (θ 0 - θ 1) * P q β s (θ 1)) / Δ q (θ 0 - θ 1)

def cramer1 (q : Model) (β s θ : Fin 2 → ℝ) : ℝ :=
  (κ q * P q β s (θ 1) - Φ q (θ 0 - θ 1) * P q β s (θ 0)) / Δ q (θ 0 - θ 1)

/-- The four scalar rows of `prop-two-student-system`:
`R_i = κ c_i + Φ(D) c_j - P(θ_i)`, `Q_0 = c_0 (A(θ_0) - c_1 H(D))`,
`Q_1 = c_1 (A(θ_1) + c_0 H(D))`, with `D = θ_0 - θ_1`. -/
def FourRows (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  (κ q * c 0 + Φ q (θ 0 - θ 1) * c 1 = P q β s (θ 0)) ∧
  (Φ q (θ 0 - θ 1) * c 0 + κ q * c 1 = P q β s (θ 1)) ∧
  (c 0 * (A q β s (θ 0) - c 1 * Hq q (θ 0 - θ 1)) = 0) ∧
  (c 1 * (A q β s (θ 1) + c 0 * Hq q (θ 0 - θ 1)) = 0)

/-- Student and teacher first moments `u = Σ c_i e(θ_i)`, `v = Σ s_k e(β_k)`
as coordinate pairs. -/
def firstMoment {n : ℕ} (c θ : Fin n → ℝ) : ℝ × ℝ :=
  (∑ i, c i * cos (θ i), ∑ i, c i * sin (θ i))

/-- First-moment matching `u = v`. -/
def FirstMomentMatched {n m : ℕ} (c θ : Fin n → ℝ) (s β : Fin m → ℝ) : Prop :=
  firstMoment c θ = firstMoment s β

/-! ## Canonical reduction of the teacher -/

/-- Projective representative in `[0, π)`. -/
def dirRep (x : ℝ) : ℝ := x - ⌊x / π⌋ * π

/-- Oriented representative in `[0, 2π)`. -/
def dirRep2 (x : ℝ) : ℝ := x - ⌊x / (2 * π)⌋ * (2 * π)

/-- The canonical centered gap `β̂ ∈ [0, π)`. -/
def canonicalGapC (β : Fin 2 → ℝ) : ℝ := dirRep (β 1 - β 0)

/-- The canonical plain-ReLU gap `β̂ ∈ [0, π]`: the oriented gap reflected into
the first half circle. -/
def canonicalGapR (β : Fin 2 → ℝ) : ℝ :=
  if dirRep2 (β 1 - β 0) ≤ π then dirRep2 (β 1 - β 0) else 2 * π - dirRep2 (β 1 - β 0)

/-- Translate every direction by `a`. -/
def translate {n : ℕ} (a : ℝ) (θ : Fin n → ℝ) : Fin n → ℝ := fun i => θ i - a

/-- Reflect every direction. -/
def reflect {n : ℕ} (θ : Fin n → ℝ) : Fin n → ℝ := fun i => -θ i

/-- Swap the two students (or teachers). -/
def swap (c : Fin 2 → ℝ) : Fin 2 → ℝ := ![c 1, c 0]

/-! ## The census -/

/-- Kind of a critical family, as read in the map legend. -/
inductive FamilyKind where
  | fit
  | colliding
  | separate
  | dead
  deriving DecidableEq

/-- Minimum type recorded by the census: the global fit, or a spurious
(nonglobal) local minimum at a split of one sign or of opposite signs. -/
inductive MinType where
  | global
  | trapPositive
  | trapMixed
  deriving DecidableEq

/-- Geometric kind of a student configuration for model `q`. -/
def familyKind (q : Model) (c θ : Fin 2 → ℝ) : FamilyKind :=
  if c 0 = 0 ∨ c 1 = 0 then .dead
  else if θ 0 ≡ θ 1 [PMOD period q] then .colliding
  else .separate

/-- The census entry of a local minimum: `(fit, global)` for a global minimum,
otherwise its kind with the sign of its split (`eq-map-census-signature`). -/
def censusEntry (q : Model) (β s c θ : Fin 2 → ℝ) : FamilyKind × MinType :=
  if IsGlobalMin q β s c θ then (.fit, .global)
  else (familyKind q c θ, if 0 < c 0 * c 1 then .trapPositive else .trapMixed)

/-- The census of a teacher: the set of entries realised by its local minima. -/
def Census (q : Model) (β s : Fin 2 → ℝ) : Set (FamilyKind × MinType) :=
  {e | ∃ c θ : Fin 2 → ℝ, IsLocalMinL q β s c θ ∧ censusEntry q β s c θ = e}

/-! ## Basic kernel calculus (proofs adapted from the paper repository) -/

theorem continuous_orthantPhi : Continuous orthantPhi := by
  unfold orthantPhi
  exact (continuous_id.mul continuous_arcsin).add
    (continuous_sqrt.comp (continuous_const.sub (continuous_id.pow 2)))

theorem orthantPhi_even (r : ℝ) : orthantPhi (-r) = orthantPhi r := by
  simp [orthantPhi]

theorem hasDerivAt_orthantPhi_of_ne {r : ℝ} (hm : r ≠ -1) (hp : r ≠ 1) :
    HasDerivAt orthantPhi (arcsin r) r := by
  have hs : 1 - r ^ 2 ≠ 0 := by
    intro h
    have : (r - 1) * (r + 1) = 0 := by nlinarith only [h]
    rcases mul_eq_zero.mp this with h | h
    · exact hp (by linarith only [h])
    · exact hm (by linarith only [h])
  have h := ((hasDerivAt_id r).mul (hasDerivAt_arcsin hm hp)).add
    ((hasDerivAt_sqrt hs).comp r ((hasDerivAt_pow 2 r).const_sub 1))
  convert h using 1; try rfl
  field_simp
  ; ring

theorem orthantPhi_of_one_le {r : ℝ} (hr : 1 ≤ r) : orthantPhi r = r * (π / 2) := by
  rw [orthantPhi, arcsin_of_one_le hr,
    sqrt_eq_zero_of_nonpos (by nlinarith only [hr] : 1 - r ^ 2 ≤ 0), add_zero]

theorem hasDerivAt_orthantPhi_one : HasDerivAt orthantPhi (π / 2) 1 := by
  have hleft : HasDerivWithinAt orthantPhi (π / 2) (Iic 1) 1 := by
    apply has_deriv_at_interval_right_endpoint_of_tendsto_deriv
      (s := Ioo (0 : ℝ) 1)
    · intro r hr
      exact (hasDerivAt_orthantPhi_of_ne (by linarith only [hr.1])
        (ne_of_lt hr.2)).differentiableAt.differentiableWithinAt
    · exact continuous_orthantPhi.continuousAt.continuousWithinAt
    · exact Ioo_mem_nhdsWithin_Iio (by norm_num : (1 : ℝ) ∈ Ioc 0 1)
    · have ht := (continuous_arcsin.tendsto (1 : ℝ)).mono_left
        (nhdsWithin_le_nhds (s := Iio 1))
      rw [arcsin_one] at ht
      apply ht.congr'
      filter_upwards [Ioo_mem_nhdsWithin_Iio
        (by norm_num : (1 : ℝ) ∈ Ioc 0 1)] with r hr
      exact (hasDerivAt_orthantPhi_of_ne (by linarith only [hr.1])
        (ne_of_lt hr.2)).deriv.symm
  have hright : HasDerivWithinAt orthantPhi (π / 2) (Ici 1) 1 := by
    have h := ((hasDerivAt_id (1 : ℝ)).mul_const (π / 2)).hasDerivWithinAt
      (s := Ici (1 : ℝ))
    simp only [one_mul] at h
    exact h.congr_of_mem (fun _ hr => orthantPhi_of_one_le hr) (by simp)
  simpa using hleft.union hright

/-- `orthantPhi' = arcsin` everywhere, including the two clamped endpoints. -/
theorem hasDerivAt_orthantPhi (r : ℝ) : HasDerivAt orthantPhi (arcsin r) r := by
  by_cases hp : r = 1
  · simpa [hp, arcsin_one] using hasDerivAt_orthantPhi_one
  by_cases hm : r = -1
  · subst r
    have hpoint : HasDerivAt orthantPhi (π / 2) (- -(1 : ℝ)) := by
      simpa using hasDerivAt_orthantPhi_one
    have h := hpoint.comp (-1) (hasDerivAt_neg' (-1))
    simpa [Function.comp_def, orthantPhi_even, arcsin_neg, arcsin_one] using h
  exact hasDerivAt_orthantPhi_of_ne hm hp

/-- `Φ_C' = -H_C` (`eq-kernel-first-derivatives`). -/
theorem hasDerivAt_phiCos (t : ℝ) : HasDerivAt phiCos (-couplingH t) t := by
  unfold phiCos couplingH
  convert (hasDerivAt_orthantPhi (cos t)).comp t (hasDerivAt_cos t) using 1
  ; ring

/-- `Φ_R' = -H_R`. -/
theorem hasDerivAt_phiCosJ (t : ℝ) : HasDerivAt phiCosJ (-couplingHJ t) t := by
  have h := (hasDerivAt_phiCos t).add ((hasDerivAt_cos t).const_mul (π / 2))
  unfold phiCosJ reluPhi couplingHJ
  convert h using 1
  ring

/-- `Φ_q' = -H_q` for both models. -/
theorem hasDerivAt_Φ (q : Model) (t : ℝ) : HasDerivAt (Φ q) (-Hq q t) t := by
  cases q
  · exact hasDerivAt_phiCos t
  · exact hasDerivAt_phiCosJ t

theorem hasDerivAt_mul_continuousAt_of_zero {f g : ℝ → ℝ} {x f' : ℝ}
    (hf : HasDerivAt f f' x) (hg : ContinuousAt g x) (hz : f x = 0) :
    HasDerivAt (fun t => f t * g t) (f' * g x) x := by
  rw [hasDerivAt_iff_tendsto_slope] at hf ⊢
  convert hf.mul (hg.mono_left nhdsWithin_le_nhds) using 1
  ext t
  simp only [slope_def_field, hz, zero_mul, sub_zero]
  ring

theorem sqrt_one_sub_cos_sq (t : ℝ) : sqrt (1 - cos t ^ 2) = |sin t| := by
  have hsq : 1 - cos t ^ 2 = sin t ^ 2 := by nlinarith only [sin_sq_add_cos_sq t]
  rw [hsq]
  exact sqrt_sq_eq_abs _

/-- The kernel equation `H_C' = Φ_C - 2|sin|`, i.e. `Φ_C'' + Φ_C = 2|sin|`
(`eq-shared-kernel-ode`), valid at every angle including the kinks. -/
theorem hasDerivAt_couplingH (t : ℝ) :
    HasDerivAt couplingH (phiCos t - 2 * |sin t|) t := by
  have key : HasDerivAt (fun x => -sin x * arcsin (cos x))
      (2 * |sin t| - phiCos t) t := by
    by_cases hz : sin t = 0
    · have h := hasDerivAt_mul_continuousAt_of_zero (hasDerivAt_sin t).neg
        (continuous_arcsin.comp continuous_cos).continuousAt (by simp [hz])
      convert h using 1
      simp only [phiCos, orthantPhi, sqrt_one_sub_cos_sq, hz, abs_zero,
        Function.comp_def]
      ring
    · have hm : cos t ≠ -1 := by
        intro h
        have : sin t ^ 2 = 0 := by nlinarith only [sin_sq_add_cos_sq t, h]
        exact hz (sq_eq_zero_iff.mp this)
      have hp : cos t ≠ 1 := by
        intro h
        have : sin t ^ 2 = 0 := by nlinarith only [sin_sq_add_cos_sq t, h]
        exact hz (sq_eq_zero_iff.mp this)
      have h := (hasDerivAt_sin t).neg.mul
        ((hasDerivAt_arcsin hm hp).comp t (hasDerivAt_cos t))
      convert h using 1
      rw [phiCos, orthantPhi, sqrt_one_sub_cos_sq]
      have habs : |sin t| ≠ 0 := (abs_pos.mpr hz).ne'
      field_simp [habs]
      nlinarith only [sq_abs (sin t)]
  have h := key.neg
  convert h using 1
  · funext x; simp only [couplingH]; ring
  · ring

/-- `H_R' = Φ_R - 2|sin|`. -/
theorem hasDerivAt_couplingHJ (t : ℝ) :
    HasDerivAt couplingHJ (phiCosJ t - 2 * |sin t|) t := by
  have h := (hasDerivAt_couplingH t).add ((hasDerivAt_sin t).const_mul (π / 2))
  unfold couplingHJ phiCosJ reluPhi
  convert h using 1
  simp only [phiCos]
  ring

/-- `H_q' = Hslope q` for both models. -/
theorem hasDerivAt_Hq (q : Model) (t : ℝ) : HasDerivAt (Hq q) (Hslope q t) t := by
  cases q
  · exact hasDerivAt_couplingH t
  · exact hasDerivAt_couplingHJ t

/-! ## One harmonic apart (`eq-planar-kernels-explicit`) -/

/-- `Φ_R = Φ_C + (π/2) cos`. -/
theorem phiCosJ_eq (x : ℝ) : phiCosJ x = phiCos x + (π / 2) * cos x := rfl

/-- `H_R = H_C + (π/2) sin`. -/
theorem couplingHJ_eq (x : ℝ) : couplingHJ x = couplingH x + (π / 2) * sin x := rfl

/-! ## Diagonal values, parity, periods -/

theorem phiCos_zero : phiCos 0 = π / 2 := by
  simp [phiCos, orthantPhi]

theorem phiCosJ_zero : phiCosJ 0 = π := by
  rw [phiCosJ_eq, phiCos_zero, cos_zero]; ring

/-- `Φ_q(0) = κ_q`. -/
theorem Φ_zero (q : Model) : Φ q 0 = κ q := by
  cases q
  · exact phiCos_zero
  · exact phiCosJ_zero

theorem phiCos_even (x : ℝ) : phiCos (-x) = phiCos x := by
  simp [phiCos, cos_neg]

theorem phiCosJ_even (x : ℝ) : phiCosJ (-x) = phiCosJ x := by
  simp [phiCosJ, reluPhi, orthantPhi, cos_neg]

theorem Φ_even (q : Model) (x : ℝ) : Φ q (-x) = Φ q x := by
  cases q
  · exact phiCos_even x
  · exact phiCosJ_even x

theorem couplingH_odd (x : ℝ) : couplingH (-x) = -couplingH x := by
  simp [couplingH, sin_neg, cos_neg]

theorem couplingHJ_odd (x : ℝ) : couplingHJ (-x) = -couplingHJ x := by
  simp [couplingHJ, couplingH_odd, sin_neg]; ring

theorem Hq_odd (q : Model) (x : ℝ) : Hq q (-x) = -Hq q x := by
  cases q
  · exact couplingH_odd x
  · exact couplingHJ_odd x

theorem phiCos_periodic : Function.Periodic phiCos π := by
  intro x
  simp [phiCos, orthantPhi, cos_add_pi]

theorem couplingH_periodic : Function.Periodic couplingH π := by
  intro x
  simp [couplingH, sin_add_pi, cos_add_pi]

theorem phiCosJ_periodic : Function.Periodic phiCosJ (2 * π) := by
  intro x
  simp [phiCosJ, reluPhi, orthantPhi, cos_add_two_pi]

theorem couplingHJ_periodic : Function.Periodic couplingHJ (2 * π) := by
  intro x
  simp [couplingHJ, couplingH, sin_add_two_pi, cos_add_two_pi]

/-- `Φ_q` has period `period q`. -/
theorem Φ_periodic (q : Model) : Function.Periodic (Φ q) (period q) := by
  cases q
  · exact phiCos_periodic
  · exact phiCosJ_periodic

theorem Hq_periodic (q : Model) : Function.Periodic (Hq q) (period q) := by
  cases q
  · exact couplingH_periodic
  · exact couplingHJ_periodic

/-- `Φ_R(π) = 0`: the plain kernel vanishes at the antipode. -/
theorem phiCosJ_pi : phiCosJ π = 0 := by
  have hπ : phiCos π = phiCos 0 := by simpa using phiCos_periodic 0
  rw [phiCosJ_eq, cos_pi, hπ, phiCos_zero]; ring

/-- `Φ_C(π/2) = 1`. -/
theorem phiCos_pi_div_two : phiCos (π / 2) = 1 := by
  simp [phiCos, orthantPhi]

/-- `H_C(π/2) = 0`. -/
theorem couplingH_pi_div_two : couplingH (π / 2) = 0 := by
  simp [couplingH]

/-- The open-branch formula `Φ_C(t) = (π/2 - t) cos t + sin t` on `[0, π]`
(`eq-centered-kernel-open-branch`). -/
theorem phiCos_branch {t : ℝ} (h0 : 0 ≤ t) (hπ : t ≤ π) :
    phiCos t = (π / 2 - t) * cos t + sin t := by
  unfold phiCos orthantPhi
  rw [sqrt_one_sub_cos_sq, abs_of_nonneg (sin_nonneg_of_nonneg_of_le_pi h0 hπ)]
  have : arcsin (cos t) = π / 2 - t := by
    rw [← sin_pi_div_two_sub, arcsin_sin] <;> linarith
  rw [this]; ring

/-- `H_C(t) = (π/2 - t) sin t` on `[0, π]`. -/
theorem couplingH_branch {t : ℝ} (h0 : 0 ≤ t) (hπ : t ≤ π) :
    couplingH t = (π / 2 - t) * sin t := by
  unfold couplingH
  have : arcsin (cos t) = π / 2 - t := by
    rw [← sin_pi_div_two_sub, arcsin_sin] <;> linarith
  rw [this]; ring

/-- `Φ_R(t) = (π - t) cos t + sin t` on `[0, π]`. -/
theorem phiCosJ_branch {t : ℝ} (h0 : 0 ≤ t) (hπ : t ≤ π) :
    phiCosJ t = (π - t) * cos t + sin t := by
  rw [phiCosJ_eq, phiCos_branch h0 hπ]; ring

/-- `H_R(t) = (π - t) sin t` on `[0, π]`. -/
theorem couplingHJ_branch {t : ℝ} (h0 : 0 ≤ t) (hπ : t ≤ π) :
    couplingHJ t = (π - t) * sin t := by
  rw [couplingHJ_eq, couplingH_branch h0 hπ]; ring

/-! ## Canonical representatives -/

theorem dirRep_nonneg (x : ℝ) : 0 ≤ dirRep x := by
  unfold dirRep
  have h := Int.floor_le (x / π)
  have : ⌊x / π⌋ * π ≤ x := by
    calc (⌊x / π⌋ : ℝ) * π ≤ (x / π) * π := by gcongr
      _ = x := by field_simp
  linarith

theorem dirRep_lt_pi (x : ℝ) : dirRep x < π := by
  unfold dirRep
  have h := Int.lt_floor_add_one (x / π)
  have : x < (⌊x / π⌋ + 1) * π := by
    calc x = (x / π) * π := by field_simp
      _ < (⌊x / π⌋ + 1) * π := by gcongr
  linarith

/-- `dirRep x ≡ x [PMOD π]`. -/
theorem dirRep_modEq (x : ℝ) : dirRep x ≡ x [PMOD π] := by
  refine ⟨⌊x / π⌋, ?_⟩
  unfold dirRep
  simp [zsmul_eq_mul]

theorem dirRep2_nonneg (x : ℝ) : 0 ≤ dirRep2 x := by
  unfold dirRep2
  have h := Int.floor_le (x / (2 * π))
  have : ⌊x / (2 * π)⌋ * (2 * π) ≤ x := by
    calc (⌊x / (2 * π)⌋ : ℝ) * (2 * π) ≤ (x / (2 * π)) * (2 * π) := by gcongr
      _ = x := by field_simp
  linarith

theorem dirRep2_lt_two_pi (x : ℝ) : dirRep2 x < 2 * π := by
  unfold dirRep2
  have h := Int.lt_floor_add_one (x / (2 * π))
  have : x < (⌊x / (2 * π)⌋ + 1) * (2 * π) := by
    calc x = (x / (2 * π)) * (2 * π) := by field_simp
      _ < (⌊x / (2 * π)⌋ + 1) * (2 * π) := by gcongr
  linarith

theorem dirRep2_modEq (x : ℝ) : dirRep2 x ≡ x [PMOD (2 * π)] := by
  refine ⟨⌊x / (2 * π)⌋, ?_⟩
  unfold dirRep2
  simp [zsmul_eq_mul]

theorem canonicalGapC_mem (β : Fin 2 → ℝ) :
    0 ≤ canonicalGapC β ∧ canonicalGapC β < π :=
  ⟨dirRep_nonneg _, dirRep_lt_pi _⟩

theorem canonicalGapR_mem (β : Fin 2 → ℝ) :
    0 ≤ canonicalGapR β ∧ canonicalGapR β ≤ π := by
  unfold canonicalGapR
  split_ifs with h
  · exact ⟨dirRep2_nonneg _, h⟩
  · push_neg at h
    exact ⟨by linarith [dirRep2_lt_two_pi (β 1 - β 0)], by linarith⟩

/-! ## Lattice helpers -/

/-- `Hq q 0 = 0`. -/
theorem Hq_zero (q : Model) : Hq q 0 = 0 := by
  cases q <;> simp [Hq, couplingH, couplingHJ]

/-- `a − b ≡ 0` iff `a ≡ b`. -/
theorem modEq_sub_zero_iff {a b p : ℝ} : a - b ≡ 0 [PMOD p] ↔ a ≡ b [PMOD p] := by
  constructor
  · rintro ⟨z, hz⟩; exact ⟨z, by linarith⟩
  · rintro ⟨z, hz⟩; exact ⟨z, by linarith⟩

/-- Collided directions see the diagonal kernel value. -/
theorem Φ_of_modEq (q : Model) {x y : ℝ} (h : x ≡ y [PMOD period q]) :
    Φ q (x - y) = κ q := by
  obtain ⟨z, hz⟩ := h
  rw [show x - y = 0 - z • period q by rw [← hz]; ring, (Φ_periodic q).sub_zsmul_eq, Φ_zero]

/-- Collided directions have zero coupling. -/
theorem Hq_of_modEq (q : Model) {x y : ℝ} (h : x ≡ y [PMOD period q]) :
    Hq q (x - y) = 0 := by
  obtain ⟨z, hz⟩ := h
  rw [show x - y = 0 - z • period q by rw [← hz]; ring, (Hq_periodic q).sub_zsmul_eq, Hq_zero]

/-- `κ_q > 0`. -/
theorem κ_pos (q : Model) : 0 < κ q := by
  cases q <;> simp [κ] <;> positivity

end TwoTeacherCensus
