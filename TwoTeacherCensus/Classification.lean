import TwoTeacherCensus.Definitions

/-!
# The two-student / two-teacher classification

Statements of the manuscript's two census theorems
(`thm:centered-two-teacher-census`, `thm:plain-two-teacher-census`; website
`thm-centered-two-student-enumeration`, `thm-plain-critical-enumeration`),
with every family of the tables made explicit, together with the proof
skeleton of Appendix C0 (Claims and Steps named in the docstrings).

Proof architecture (see `PROOF_PLAN.md`):
1. **Reduction** (`stp:census-centered-reduction`, `stp:census-plain-reduction`):
   exact isometries of the loss and the canonical gap.
2. **Rows** (`prop-two-student-system`): criticality is the four scalar rows.
3. **First moment** (`eq-general-harmonic-split`): `L_R = L_C + ⅛‖u-v‖²`.
4. **Collision block** (`prop-collision-block`, `clm:census-centered-selector`,
   `lem-flat-sufficiency`): the split line is flat, rotation reads the teacher
   potential, and the angular Hessian selector decides the type.
5. **Case split** by strata (zero field / one line / genuine) and by student
   precedence (both dead / one dead / collided / separated).
6. **Second variation** on the separated stratum (beam, quarter, Schur).

Lemmas that are short algebra are proved here; everything else carries a
`sorry` and a docstring naming the Claim/Step it realises.
-/

noncomputable section

open Real Set Filter Topology Classical
open scoped BigOperators

namespace TwoTeacherCensus

/-! ## 1. Exact isometries of the loss (Step `stp:census-centered-reduction`,
`stp:census-plain-reduction`) -/

section Symmetries

variable {n m : ℕ} (q : Model)

/-- Translating every direction leaves the loss unchanged (`Φ_q` depends on
differences only). -/
theorem L_translate (a : ℝ) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    L q (translate a β) s c (translate a θ) = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy translate
  simp only [sub_sub_sub_cancel_right]

/-- Reflecting every direction leaves the loss unchanged (`Φ_q` is even). -/
theorem L_reflect (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    L q (reflect β) s c (reflect θ) = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy reflect
  have h : ∀ x y : ℝ, Φ q (-x - -y) = Φ q (x - y) := fun x y => by
    rw [show -x - -y = -(x - y) by ring, Φ_even]
  simp only [h]

/-- Negating all masses leaves the loss unchanged (it is a quadratic form in
the concatenated mass vector). -/
theorem L_neg_masses (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    L q β (-s) (-c) θ = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy
  simp only [Pi.neg_apply, neg_mul_neg]

/-- The two students may be relabelled. -/
theorem L_swap_students (β s : Fin m → ℝ) (c θ : Fin 2 → ℝ) :
    L q β s (swap c) (swap θ) = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy swap
  have h : Φ q (θ 1 - θ 0) = Φ q (θ 0 - θ 1) := by
    rw [← Φ_even, neg_sub]
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.head_cons, h]
  ring

/-- The two teachers may be relabelled. -/
theorem L_swap_teachers (β s : Fin 2 → ℝ) (c θ : Fin n → ℝ) :
    L q (swap β) (swap s) c θ = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy swap
  have h : Φ q (β 1 - β 0) = Φ q (β 0 - β 1) := by
    rw [← Φ_even, neg_sub]
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.head_cons, h, Finset.sum_add_distrib]
  ring

theorem Φ_update_period (θ : Fin n → ℝ) (i0 i j : Fin n) :
    Φ q (Function.update θ i0 (θ i0 + period q) i -
        Function.update θ i0 (θ i0 + period q) j) = Φ q (θ i - θ j) := by
  rcases eq_or_ne i i0 with hi | hi <;> rcases eq_or_ne j i0 with hj | hj
  · rw [hi, hj, Function.update_same]; simp
  · rw [hi, Function.update_same, Function.update_noteq hj,
      show θ i0 + period q - θ j = (θ i0 - θ j) + period q by ring, Φ_periodic q]
  · rw [hj, Function.update_same, Function.update_noteq hi,
      show θ i - (θ i0 + period q) = (θ i - θ i0) - period q by ring, (Φ_periodic q).sub_eq]
  · rw [Function.update_noteq hi, Function.update_noteq hj]

theorem Φ_update_period_teacher (θ : Fin n → ℝ) (β : Fin m → ℝ) (i0 i : Fin n) (k : Fin m) :
    Φ q (Function.update θ i0 (θ i0 + period q) i - β k) = Φ q (θ i - β k) := by
  rcases eq_or_ne i i0 with hi | hi
  · rw [hi, Function.update_same,
      show θ i0 + period q - β k = (θ i0 - β k) + period q by ring, Φ_periodic q]
  · rw [Function.update_noteq hi]

/-- Shifting one student direction by the kernel period (`π` for centered, `2π`
for plain ReLU) leaves the loss unchanged.  For the centered model this is the
orientation-fiber symmetry of `one-harmonic-apart`. -/
theorem L_shift_period (β s : Fin m → ℝ) (c θ : Fin n → ℝ) (i0 : Fin n) :
    L q β s c (Function.update θ i0 (θ i0 + period q)) = L q β s c θ := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy
  simp only [Φ_update_period]
  simp only [Φ_update_period_teacher]

/-- The Fréchet derivative of a translate of the argument. -/
theorem fderiv_comp_add_const' {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (g : E → ℝ) (v x : E) :
    fderiv ℝ (fun p => g (p + v)) x = fderiv ℝ g (x + v) := by
  by_cases hd : DifferentiableAt ℝ g (x + v)
  · have h1 : HasFDerivAt (fun p : E => p + v) (ContinuousLinearMap.id ℝ E) x :=
      (hasFDerivAt_id x).add_const v
    have h2 := hd.hasFDerivAt.comp x h1
    rw [ContinuousLinearMap.comp_id] at h2
    exact h2.fderiv
  · have hd' : ¬ DifferentiableAt ℝ (fun p => g (p + v)) x := by
      intro h
      apply hd
      have h2 : DifferentiableAt ℝ (fun p : E => p - v) (x + v) :=
        (differentiableAt_id.sub_const v)
      have h3 : DifferentiableAt ℝ (fun p => g (p + v)) ((x + v) - v) := by
        simpa using h
      have := h3.comp (x + v) h2
      convert this using 1
      funext p
      simp
    rw [fderiv_zero_of_not_differentiableAt hd, fderiv_zero_of_not_differentiableAt hd']

/-- Local minimality of a translate of the argument. -/
theorem isLocalMin_comp_add_const' {E : Type*} [NormedAddCommGroup E] [NormedSpace ℝ E]
    (g : E → ℝ) (v x : E) :
    IsLocalMin (fun p => g (p + v)) x ↔ IsLocalMin g (x + v) := by
  constructor
  · intro h
    have h0 : IsLocalMin (fun p => g (p + v)) ((x + v) - v) := by simpa using h
    have h' := h0.comp_continuous (g := fun p : E => p - v) (b := x + v)
      (continuous_id.sub continuous_const).continuousAt
    convert h' using 1
    funext p
    simp
  · intro h
    exact h.comp_continuous (g := fun p : E => p + v) (b := x)
      (continuous_id.add continuous_const).continuousAt

/-- The translated loss is the loss at a translated argument. -/
theorem lossPair_translate_eq (a : ℝ) (β s : Fin m → ℝ) :
    lossPair q (translate a β) s =
      fun p : (Fin n → ℝ) × (Fin n → ℝ) =>
        lossPair q β s (p + ((0 : Fin n → ℝ), fun _ => a)) := by
  funext p
  unfold lossPair
  have hθ : translate a (fun i => p.2 i + a) = p.2 := by
    funext i; simp [translate]
  have h := L_translate q a β s p.1 (fun i => p.2 i + a)
  rw [hθ] at h
  rw [h]
  simp only [Prod.fst_add, Prod.snd_add, add_zero]
  rfl

/-- Criticality transports along translation.  Realises the translation clause
of `stp:census-centered-reduction` (long tree: `critical_translate_iff`). -/
theorem isCritical_translate_iff (a : ℝ) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    IsCritical q (translate a β) s c (translate a θ) ↔ IsCritical q β s c θ := by
  unfold IsCritical
  rw [lossPair_translate_eq, fderiv_comp_add_const']
  have : ((c, translate a θ) + ((0 : Fin n → ℝ), fun _ => a)) = (c, θ) := by
    ext i <;> simp [translate]
  rw [this]

/-- Local minimality transports along translation. -/
theorem isLocalMinL_translate_iff (a : ℝ) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    IsLocalMinL q (translate a β) s c (translate a θ) ↔ IsLocalMinL q β s c θ := by
  unfold IsLocalMinL
  rw [lossPair_translate_eq, isLocalMin_comp_add_const']
  have : ((c, translate a θ) + ((0 : Fin n → ℝ), fun _ => a)) = (c, θ) := by
    ext i <;> simp [translate]
  rw [this]

/-- Global minimality transports along translation (proved: the loss values
agree pointwise and translation is a bijection of the student angles). -/
theorem isGlobalMin_translate_iff (a : ℝ) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    IsGlobalMin q (translate a β) s c (translate a θ) ↔ IsGlobalMin q β s c θ := by
  unfold IsGlobalMin
  constructor
  · intro h c' θ'
    have := h c' (translate a θ')
    rwa [L_translate, L_translate] at this
  · intro h c' θ'
    rw [L_translate]
    have hθ : θ' = translate a (translate (-a) θ') := by
      funext i; simp [translate]
    calc L q β s c θ ≤ L q β s c' (translate (-a) θ') := h _ _
      _ = L q (translate a β) s c' (translate a (translate (-a) θ')) :=
          (L_translate q a β s c' _).symm
      _ = L q (translate a β) s c' θ' := by rw [← hθ]

end Symmetries

/-! ## 2. The teacher strata -/

section Strata

/-- Centered zero field: both masses zero, or two coincident lines with
cancelling masses. -/
def CenteredZeroField (β s : Fin 2 → ℝ) : Prop :=
  (s 0 = 0 ∧ s 1 = 0) ∨ (β 0 ≡ β 1 [PMOD π] ∧ s 0 + s 1 = 0)

/-- Centered teacher with exactly one effective line at `a` of mass `σ ≠ 0`:
coincident lines with `σ = s_0 + s_1`, or exactly one raw mass surviving. -/
def CenteredOneLine (β s : Fin 2 → ℝ) (a σ : ℝ) : Prop :=
  σ ≠ 0 ∧
    ((β 0 ≡ β 1 [PMOD π] ∧ a = β 0 ∧ σ = s 0 + s 1) ∨
      (s 1 = 0 ∧ a = β 0 ∧ σ = s 0) ∨
      (s 0 = 0 ∧ a = β 1 ∧ σ = s 1))

/-- Genuine centered teacher: two distinct lines, both masses nonzero. -/
def CenteredGenuine (β s : Fin 2 → ℝ) : Prop :=
  ¬ β 0 ≡ β 1 [PMOD π] ∧ s 0 ≠ 0 ∧ s 1 ≠ 0

/-- The centered strata are exhaustive. -/
theorem centered_strata_exhaustive (β s : Fin 2 → ℝ) :
    CenteredZeroField β s ∨ (∃ a σ, CenteredOneLine β s a σ) ∨ CenteredGenuine β s := by
  by_cases hβ : β 0 ≡ β 1 [PMOD π]
  · by_cases hsum : s 0 + s 1 = 0
    · exact Or.inl (Or.inr ⟨hβ, hsum⟩)
    · exact Or.inr (Or.inl ⟨β 0, s 0 + s 1, hsum, Or.inl ⟨hβ, rfl, rfl⟩⟩)
  · by_cases h0 : s 0 = 0
    · by_cases h1 : s 1 = 0
      · exact Or.inl (Or.inl ⟨h0, h1⟩)
      · exact Or.inr (Or.inl ⟨β 1, s 1, h1, Or.inr (Or.inr ⟨h0, rfl, rfl⟩)⟩)
    · by_cases h1 : s 1 = 0
      · exact Or.inr (Or.inl ⟨β 0, s 0, h0, Or.inr (Or.inl ⟨h1, rfl, rfl⟩)⟩)
      · exact Or.inr (Or.inr ⟨hβ, h0, h1⟩)

/-- Plain-ReLU zero field: both masses zero, or coincident rays with
cancelling masses. -/
def PlainZeroField (β s : Fin 2 → ℝ) : Prop :=
  (s 0 = 0 ∧ s 1 = 0) ∨ (β 0 ≡ β 1 [PMOD (2 * π)] ∧ s 0 + s 1 = 0)

/-- Plain-ReLU teacher with one effective ray at `a` of mass `σ ≠ 0`. -/
def PlainOneRay (β s : Fin 2 → ℝ) (a σ : ℝ) : Prop :=
  σ ≠ 0 ∧
    ((β 0 ≡ β 1 [PMOD (2 * π)] ∧ a = β 0 ∧ σ = s 0 + s 1) ∨
      (s 1 = 0 ∧ a = β 0 ∧ σ = s 0) ∨
      (s 0 = 0 ∧ a = β 1 ∧ σ = s 1))

/-- Antipodal plain-ReLU teacher: `β̂ = π`, both masses nonzero. -/
def PlainAntipodal (β s : Fin 2 → ℝ) : Prop :=
  β 1 ≡ β 0 + π [PMOD (2 * π)] ∧ s 0 ≠ 0 ∧ s 1 ≠ 0

/-- Genuine plain-ReLU teacher: `0 < β̂ < π`, both masses nonzero. -/
def PlainGenuine (β s : Fin 2 → ℝ) : Prop :=
  ¬ β 1 ≡ β 0 [PMOD (2 * π)] ∧ ¬ β 1 ≡ β 0 + π [PMOD (2 * π)] ∧ s 0 ≠ 0 ∧ s 1 ≠ 0

/-- Two distinct rays with nonzero masses: antipodal or genuine. -/
def PlainTwoRay (β s : Fin 2 → ℝ) : Prop := PlainAntipodal β s ∨ PlainGenuine β s

/-- The plain-ReLU strata are exhaustive. -/
theorem plain_strata_exhaustive (β s : Fin 2 → ℝ) :
    PlainZeroField β s ∨ (∃ a σ, PlainOneRay β s a σ) ∨ PlainAntipodal β s ∨ PlainGenuine β s := by
  by_cases hβ : β 0 ≡ β 1 [PMOD (2 * π)]
  · by_cases hsum : s 0 + s 1 = 0
    · exact Or.inl (Or.inr ⟨hβ, hsum⟩)
    · exact Or.inr (Or.inl ⟨β 0, s 0 + s 1, hsum, Or.inl ⟨hβ, rfl, rfl⟩⟩)
  · by_cases h0 : s 0 = 0
    · by_cases h1 : s 1 = 0
      · exact Or.inl (Or.inl ⟨h0, h1⟩)
      · exact Or.inr (Or.inl ⟨β 1, s 1, h1, Or.inr (Or.inr ⟨h0, rfl, rfl⟩)⟩)
    · by_cases h1 : s 1 = 0
      · exact Or.inr (Or.inl ⟨β 0, s 0, h0, Or.inr (Or.inl ⟨h1, rfl, rfl⟩)⟩)
      · by_cases hanti : β 1 ≡ β 0 + π [PMOD (2 * π)]
        · exact Or.inr (Or.inr (Or.inl ⟨hanti, h0, h1⟩))
        · exact Or.inr (Or.inr (Or.inr ⟨fun h => hβ h.symm, hanti, h0, h1⟩))

/-- A centered zero field has zero teacher potential. -/
theorem P_eq_zero_of_centeredZeroField {β s : Fin 2 → ℝ} (h : CenteredZeroField β s) (t : ℝ) :
    P .centered β s t = 0 := by
  unfold P
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨⟨z, hz⟩, hsum⟩
  · simp [h0, h1]
  · have hΦ : Φ .centered (t - β 1) = Φ .centered (t - β 0) := by
      rw [show t - β 1 = (t - β 0) - z • period .centered by
        simp only [period]; rw [← hz]; ring]
      exact (Φ_periodic .centered).sub_zsmul_eq z
    rw [hΦ, show s 1 = -s 0 by linarith]
    ring

/-- A centered zero field has zero torque. -/
theorem A_eq_zero_of_centeredZeroField {β s : Fin 2 → ℝ} (h : CenteredZeroField β s) (t : ℝ) :
    A .centered β s t = 0 := by
  unfold A
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨⟨z, hz⟩, hsum⟩
  · simp [h0, h1]
  · have hH : Hq .centered (t - β 1) = Hq .centered (t - β 0) := by
      rw [show t - β 1 = (t - β 0) - z • period .centered by
        simp only [period]; rw [← hz]; ring]
      exact (Hq_periodic .centered).sub_zsmul_eq z
    rw [hH, show s 1 = -s 0 by linarith]
    ring

/-- A centered zero field has zero self-energy. -/
theorem selfEnergy_eq_zero_of_centeredZeroField {β s : Fin 2 → ℝ} (h : CenteredZeroField β s) :
    kernelTeacherSelfEnergy (Φ .centered) β s = 0 := by
  unfold kernelTeacherSelfEnergy
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨hβ, hsum⟩
  · simp [h0, h1]
  · have hev : Φ .centered (β 1 - β 0) = Φ .centered (β 0 - β 1) := by rw [← Φ_even, neg_sub]
    rw [hev, Φ_of_modEq .centered hβ, sub_self, Φ_zero, sub_self, Φ_zero,
      show s 1 = -s 0 by linarith]
    ring

/-- A plain zero field has zero self-energy. -/
theorem selfEnergy_eq_zero_of_plainZeroField {β s : Fin 2 → ℝ} (h : PlainZeroField β s) :
    kernelTeacherSelfEnergy (Φ .plainRelu) β s = 0 := by
  unfold kernelTeacherSelfEnergy
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨hβ, hsum⟩
  · simp [h0, h1]
  · have hev : Φ .plainRelu (β 1 - β 0) = Φ .plainRelu (β 0 - β 1) := by rw [← Φ_even, neg_sub]
    rw [hev, Φ_of_modEq .plainRelu hβ, sub_self, Φ_zero, sub_self, Φ_zero,
      show s 1 = -s 0 by linarith]
    ring

/-- A plain zero field has zero teacher potential. -/
theorem P_eq_zero_of_plainZeroField {β s : Fin 2 → ℝ} (h : PlainZeroField β s) (t : ℝ) :
    P .plainRelu β s t = 0 := by
  unfold P
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨⟨z, hz⟩, hsum⟩
  · simp [h0, h1]
  · have hΦ : Φ .plainRelu (t - β 1) = Φ .plainRelu (t - β 0) := by
      rw [show t - β 1 = (t - β 0) - z • period .plainRelu by
        simp only [period]; rw [← hz]; ring]
      exact (Φ_periodic .plainRelu).sub_zsmul_eq z
    rw [hΦ, show s 1 = -s 0 by linarith]
    ring

/-- A plain zero field has zero torque. -/
theorem A_eq_zero_of_plainZeroField {β s : Fin 2 → ℝ} (h : PlainZeroField β s) (t : ℝ) :
    A .plainRelu β s t = 0 := by
  unfold A
  simp only [Fin.sum_univ_two]
  rcases h with ⟨h0, h1⟩ | ⟨⟨z, hz⟩, hsum⟩
  · simp [h0, h1]
  · have hH : Hq .plainRelu (t - β 1) = Hq .plainRelu (t - β 0) := by
      rw [show t - β 1 = (t - β 0) - z • period .plainRelu by
        simp only [period]; rw [← hz]; ring]
      exact (Hq_periodic .plainRelu).sub_zsmul_eq z
    rw [hH, show s 1 = -s 0 by linarith]
    ring

end Strata

/-! ## 3. One harmonic apart: `L_R = L_C + ⅛‖u − v‖²` (`eq-general-harmonic-split`) -/

section FirstMoment

theorem double_sum_add_cos {N M : ℕ} (kf : ℝ → ℝ) (r : ℝ) (w φ : Fin N → ℝ)
    (w' φ' : Fin M → ℝ) :
    ∑ i, ∑ j, w i * w' j * (kf (φ i - φ' j) + r * cos (φ i - φ' j)) =
      ∑ i, ∑ j, w i * w' j * kf (φ i - φ' j) +
        r * ((∑ i, w i * cos (φ i)) * (∑ j, w' j * cos (φ' j)) +
          (∑ i, w i * sin (φ i)) * (∑ j, w' j * sin (φ' j))) := by
  have e1 : (∑ i, w i * cos (φ i)) * (∑ j, w' j * cos (φ' j)) =
      ∑ i, ∑ j, w i * cos (φ i) * (w' j * cos (φ' j)) := by
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl fun i _ => by rw [Finset.mul_sum]
  have e2 : (∑ i, w i * sin (φ i)) * (∑ j, w' j * sin (φ' j)) =
      ∑ i, ∑ j, w i * sin (φ i) * (w' j * sin (φ' j)) := by
    rw [Finset.sum_mul]
    exact Finset.sum_congr rfl fun i _ => by rw [Finset.mul_sum]
  rw [e1, e2, ← Finset.sum_add_distrib, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [← Finset.sum_add_distrib, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [cos_sub]; ring

/-- Adding a first harmonic `r cos` to a kernel adds `(r/2)‖u − v‖²` to the
excess loss, `u, v` the student and teacher first moments. -/
theorem kernelExcessLoss_add_cos {n m : ℕ} (kf : ℝ → ℝ) (r : ℝ) (c θ : Fin n → ℝ)
    (β s : Fin m → ℝ) :
    kernelExcessLoss (fun x => kf x + r * cos x) c β s θ =
      kernelExcessLoss kf c β s θ +
        (r / 2) * (((firstMoment c θ).1 - (firstMoment s β).1) ^ 2 +
          ((firstMoment c θ).2 - (firstMoment s β).2) ^ 2) := by
  unfold kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy firstMoment
  simp only [double_sum_add_cos]
  ring

/-- **The mean-residual decomposition** `L_R = L_C + ⅛‖u − v‖²`
(`eq-general-harmonic-split`; long tree `L_R_eq_L_C_add_firstMomentMismatchSq`),
for arbitrary widths. -/
theorem L_R_eq_L_C_add_firstMoment {n m : ℕ} (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    L_R β s c θ = L_C β s c θ +
      (1 / 8) * (((firstMoment c θ).1 - (firstMoment s β).1) ^ 2 +
        ((firstMoment c θ).2 - (firstMoment s β).2) ^ 2) := by
  have hΦ : Φ .plainRelu = fun x => Φ .centered x + (π / 2) * cos x := by
    funext x; exact phiCosJ_eq x
  unfold L_R L_C L
  rw [hΦ, kernelExcessLoss_add_cos]
  field_simp
  ring

/-- On the live-student locus, criticality for both losses is criticality for
one of them plus first-moment matching (long tree
`critical_both_iff_centered_and_matched`). -/
theorem isCritical_both_iff_matched {β s c θ : Fin 2 → ℝ} (halive : ∃ i, c i ≠ 0) :
    (IsCritical .centered β s c θ ∧ IsCritical .plainRelu β s c θ) ↔
      (IsCritical .centered β s c θ ∧ FirstMomentMatched c θ s β) := by
  sorry

end FirstMoment

/-! ## 4. Nonnegativity via the Gaussian representation -/

section Gaussian

/-- The planar standard Gaussian density. -/
def stdGaussianDensity2 (x : EuclideanSpace ℝ (Fin 2)) : ℝ :=
  (sqrt (2 * π))⁻¹ ^ 2 * exp (-‖x‖ ^ 2 / 2)

/-- The unit direction `e(θ) = (cos θ, sin θ)`. -/
def Angle (θ : ℝ) : EuclideanSpace ℝ (Fin 2) := ![cos θ, sin θ]

/-- The two features `C(z) = |z|/2` and `R(z) = max 0 z`. -/
def Feature : Model → ℝ → ℝ
  | .centered, z => |z| / 2
  | .plainRelu, z => max 0 z

/-- The literal planar population loss `½ E[(student − teacher)²]`. -/
def GaussianLoss {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) : ℝ :=
  (1 / 2 : ℝ) * ∫ x : EuclideanSpace ℝ (Fin 2),
    ((∑ i, c i * Feature q (inner (Angle (θ i)) x)) -
      ∑ k, s k * Feature q (inner (Angle (β k)) x)) ^ 2 * stdGaussianDensity2 x

theorem stdGaussianDensity2_nonneg (x : EuclideanSpace ℝ (Fin 2)) :
    0 ≤ stdGaussianDensity2 x := by
  unfold stdGaussianDensity2
  positivity

theorem gaussianLoss_nonneg {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    0 ≤ GaussianLoss q β s c θ := by
  unfold GaussianLoss
  apply mul_nonneg (by norm_num)
  apply MeasureTheory.integral_nonneg
  intro x
  exact mul_nonneg (sq_nonneg _) (stdGaussianDensity2_nonneg x)

/-- **The kernel form of the population loss** (`eq-loss-quadratic`,
`lem:pair-moments`): the finite kernel energy `L` is the literal Gaussian
loss.  This is the only place where Gaussian integration enters; the paper
repository proves it as `eq_loss_quadratic` via the pair moment
`E|⟨u,x⟩||⟨z,x⟩| = (2/π)(ρ arcsin ρ + √(1-ρ²))`.  Port of ≈ 1300 lines of
`Preliminaries.lean` (rotation invariance, polar pair moment). -/
theorem L_eq_gaussianLoss {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    L q β s c θ = GaussianLoss q β s c θ := by
  sorry

/-- The population loss is nonnegative. -/
theorem L_nonneg {n m : ℕ} (q : Model) (β s : Fin m → ℝ) (c θ : Fin n → ℝ) :
    0 ≤ L q β s c θ := by
  rw [L_eq_gaussianLoss]; exact gaussianLoss_nonneg q β s c θ

/-- Zero loss is global minimality. -/
theorem isGlobalMin_of_L_eq_zero {n m : ℕ} {q : Model} {β s : Fin m → ℝ} {c θ : Fin n → ℝ}
    (h : L q β s c θ = 0) : IsGlobalMin q β s c θ := by
  intro c' θ'
  rw [h]; exact L_nonneg q β s c' θ'

end Gaussian

/-! ## 5. Four rows are exactly stationarity (`prop-two-student-system`) -/

section Rows

variable {β s c θ : Fin 2 → ℝ} {q : Model}

/-- The loss is differentiable in the joint parameter (the kernels are `C¹`). -/
theorem differentiableAt_lossPair (q : Model) (β s : Fin 2 → ℝ) (p : (Fin 2 → ℝ) × (Fin 2 → ℝ)) :
    DifferentiableAt ℝ (lossPair q β s) p := by
  sorry

/-- **Criticality is exactly the four scalar rows** (`prop-two-student-system`;
long tree `isCriticalL_C_two_student_iff_four_equations`,
`isCriticalL_R_two_student_iff_four_equations`):
`2π ∂L/∂c_i = R_i` and `2π ∂L/∂θ_i = Q_i`. -/
theorem isCritical_iff_fourRows (q : Model) (β s c θ : Fin 2 → ℝ) :
    IsCritical q β s c θ ↔ FourRows q β s c θ := by
  sorry


/-- The exact mass curve (`eq-classification-mass-curve`): moving one mass by
`r` changes the loss by `(κ_q/(4π)) r² + (r/(2π)) R_0`, exactly. -/
theorem L_mass0_curve (q : Model) (β s c θ : Fin 2 → ℝ) (r : ℝ) :
    L q β s (Function.update c 0 (c 0 + r)) θ - L q β s c θ =
      (κ q / (4 * π)) * r ^ 2 +
        (r / (2 * π)) * (κ q * c 0 + Φ q (θ 0 - θ 1) * c 1 - P q β s (θ 0)) := by
  have hev : Φ q (θ 1 - θ 0) = Φ q (θ 0 - θ 1) := by rw [← Φ_even, neg_sub]
  have h0 : Φ q (θ 0 - θ 0) = κ q := by rw [sub_self, Φ_zero]
  have h1 : Φ q (θ 1 - θ 1) = κ q := by rw [sub_self, Φ_zero]
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy P
  simp only [Fin.sum_univ_two, Function.update_same,
    Function.update_noteq (show (1 : Fin 2) ≠ 0 by decide), hev, h0, h1]
  field_simp
  ring

/-- **No critical point is a local maximum** (`lem-no-local-max`): along the
mass curve the loss increases quadratically. -/
theorem not_isLocalMaxL_of_fourRows (h : FourRows q β s c θ) :
    ¬ IsLocalMaxL q β s c θ := by
  intro hmax
  have hR0 : κ q * c 0 + Φ q (θ 0 - θ 1) * c 1 - P q β s (θ 0) = 0 := by linarith [h.1]
  have hcont : Continuous fun r : ℝ =>
      ((Function.update c 0 (c 0 + r), θ) : (Fin 2 → ℝ) × (Fin 2 → ℝ)) := by
    apply Continuous.prod_mk _ continuous_const
    exact Continuous.update continuous_const 0 (continuous_const.add continuous_id)
  have h0 : ((Function.update c 0 (c 0 + 0), θ) : (Fin 2 → ℝ) × (Fin 2 → ℝ)) = (c, θ) := by
    simp
  have hmax' : IsLocalMax
      (lossPair q β s ∘ fun r : ℝ => ((Function.update c 0 (c 0 + r), θ))) 0 := by
    apply IsLocalMax.comp_continuous _ hcont.continuousAt
    rw [h0]; exact hmax
  rcases Metric.eventually_nhds_iff.mp hmax' with ⟨ε, hε, hball⟩
  have hr : dist (ε / 2) (0 : ℝ) < ε := by
    rw [Real.dist_eq, sub_zero, abs_of_pos (by positivity)]; linarith
  have hle := hball hr
  simp only [Function.comp, lossPair, add_zero, Function.update_eq_self] at hle
  have hcurve := L_mass0_curve q β s c θ (ε / 2)
  rw [hR0, mul_zero, add_zero] at hcurve
  have hκ : 0 < κ q := by cases q <;> simp [κ] <;> positivity
  have hpos : 0 < κ q / (4 * π) * (ε / 2) ^ 2 := by positivity
  linarith

theorem not_isLocalMaxL_of_isCritical (h : IsCritical q β s c θ) :
    ¬ IsLocalMaxL q β s c θ :=
  not_isLocalMaxL_of_fourRows ((isCritical_iff_fourRows q β s c θ).mp h)

/-- At a critical point, "topological saddle" is exactly "not a local minimum". -/
theorem isTopologicalSaddle_iff_not_isLocalMin (h : IsCritical q β s c θ) :
    IsTopologicalSaddle q β s c θ ↔ ¬ IsLocalMinL q β s c θ := by
  constructor
  · exact fun hs => hs.2.1
  · exact fun hmin => ⟨h, hmin, not_isLocalMaxL_of_isCritical h⟩

/-- Every critical point carries exactly one of the three labels. -/
theorem variationalType_of_isCritical (h : IsCritical q β s c θ) :
    VariationalType q β s c θ := by
  unfold VariationalType
  by_cases hmin : IsLocalMinL q β s c θ
  · by_cases hglob : IsGlobalMin q β s c θ
    · exact Or.inl hglob
    · exact Or.inr (Or.inl ⟨hmin, hglob⟩)
  · exact Or.inr (Or.inr ((isTopologicalSaddle_iff_not_isLocalMin h).mpr hmin))

/-- A global minimum is a local minimum. -/
theorem isLocalMinL_of_isGlobalMin {n m : ℕ} {q : Model} {β s : Fin m → ℝ} {c θ : Fin n → ℝ}
    (h : IsGlobalMin q β s c θ) : IsLocalMinL q β s c θ :=
  Filter.eventually_of_forall fun p => h p.1 p.2

/-- The three labels are pairwise exclusive. -/
theorem variationalType_exclusive :
    ¬ (IsGlobalMin q β s c θ ∧ IsSpuriousMin q β s c θ) ∧
    ¬ (IsGlobalMin q β s c θ ∧ IsTopologicalSaddle q β s c θ) ∧
    ¬ (IsSpuriousMin q β s c θ ∧ IsTopologicalSaddle q β s c θ) := by
  refine ⟨fun h => h.2.2 h.1, fun h => h.2.2.1 (isLocalMinL_of_isGlobalMin h.1),
    fun h => h.2.2.1 h.1.1⟩

/-- The kernel window `|Φ_q(D)| < κ_q` off the lattice `period q · ℤ`
(`prop-kernel-properties`(6)), hence `Δ_q(D) > 0`. -/
theorem Δ_pos_of_not_modEq (q : Model) {D : ℝ} (hD : ¬ D ≡ 0 [PMOD period q]) :
    0 < Δ q D := by
  sorry

/-- `Δ_q(D) = 0` exactly on the lattice (so a separated pair always has an
invertible radial block). -/
theorem Δ_eq_zero_iff (q : Model) (D : ℝ) : Δ q D = 0 ↔ D ≡ 0 [PMOD period q] := by
  sorry

/-- Cramér's rule (`eq-separated-cramer`): with `Δ_q(D) ≠ 0` the two radial
rows are equivalent to the two displayed mass formulas. -/
theorem radialRows_iff_cramer (hΔ : Δ q (θ 0 - θ 1) ≠ 0) :
    ((κ q * c 0 + Φ q (θ 0 - θ 1) * c 1 = P q β s (θ 0)) ∧
      (Φ q (θ 0 - θ 1) * c 0 + κ q * c 1 = P q β s (θ 1))) ↔
    (c 0 = cramer0 q β s θ ∧ c 1 = cramer1 q β s θ) := by
  unfold cramer0 cramer1
  constructor
  · rintro ⟨h0, h1⟩
    constructor
    · rw [eq_div_iff hΔ]; unfold Δ at hΔ ⊢; linear_combination κ q * h0 - Φ q (θ 0 - θ 1) * h1
    · rw [eq_div_iff hΔ]; unfold Δ at hΔ ⊢; linear_combination κ q * h1 - Φ q (θ 0 - θ 1) * h0
  · rintro ⟨h0, h1⟩
    rw [h0, h1]
    unfold Δ at hΔ ⊢
    constructor
    · field_simp; ring
    · field_simp; ring

/-- The angular rows solve the masses a second time when both students are
live and `H_q(D) ≠ 0` (`eq-classification-angular-derivative-masses`). -/
theorem angularRows_iff_torqueMasses (hc0 : c 0 ≠ 0) (hc1 : c 1 ≠ 0)
    (hH : Hq q (θ 0 - θ 1) ≠ 0) :
    ((c 0 * (A q β s (θ 0) - c 1 * Hq q (θ 0 - θ 1)) = 0) ∧
      (c 1 * (A q β s (θ 1) + c 0 * Hq q (θ 0 - θ 1)) = 0)) ↔
    (c 1 = A q β s (θ 0) / Hq q (θ 0 - θ 1) ∧
      c 0 = - A q β s (θ 1) / Hq q (θ 0 - θ 1)) := by
  constructor
  · rintro ⟨h0, h1⟩
    rcases mul_eq_zero.mp h0 with h | h
    · exact absurd h hc0
    rcases mul_eq_zero.mp h1 with h' | h'
    · exact absurd h' hc1
    constructor
    · rw [eq_div_iff hH]; linarith
    · rw [eq_div_iff hH]; linarith
  · rintro ⟨h1, h0⟩
    constructor
    · rw [h1]; field_simp
    · rw [h0]; field_simp

/-- A zero teacher field (`P_q ≡ 0 ≡ A_q`) has exactly the zero student
measures as critical points: `c = 0`, or a collided cancelling pair.
Realises the zero-field rows of `stp:census-centered-one-line` /
`stp:census-plain-reduction` for both models at once. -/
theorem zeroField_critical_iff (hP : ∀ t, P q β s t = 0) (hA : ∀ t, A q β s t = 0) :
    IsCritical q β s c θ ↔
      ((c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD period q] ∧ c 0 + c 1 = 0)) := by
  rw [isCritical_iff_fourRows]
  unfold FourRows
  rw [hP, hP, hA, hA]
  by_cases hD : θ 0 ≡ θ 1 [PMOD period q]
  · rw [Φ_of_modEq q hD, Hq_of_modEq q hD]
    have hκ := (κ_pos q).ne'
    constructor
    · rintro ⟨_, h1, _, _⟩
      right
      refine ⟨hD, ?_⟩
      have : κ q * (c 0 + c 1) = 0 := by linarith
      exact (mul_eq_zero.mp this).resolve_left hκ
    · rintro (⟨h0, h1⟩ | ⟨_, hsum⟩)
      · simp [h0, h1]
      · refine ⟨by linear_combination κ q * hsum, by linear_combination κ q * hsum, ?_, ?_⟩ <;>
          simp
  · have hΔ : Δ q (θ 0 - θ 1) ≠ 0 :=
      (Δ_pos_of_not_modEq q (fun h' => hD (modEq_sub_zero_iff.mp h'))).ne'
    constructor
    · rintro ⟨h0, h1, _, _⟩
      left
      have h0' : κ q * c 0 + Φ q (θ 0 - θ 1) * c 1 = P q β s (θ 0) := by rw [hP]; exact h0
      have h1' : Φ q (θ 0 - θ 1) * c 0 + κ q * c 1 = P q β s (θ 1) := by rw [hP]; exact h1
      have hc := (radialRows_iff_cramer (q := q) (β := β) (s := s) hΔ).mp ⟨h0', h1'⟩
      unfold cramer0 cramer1 at hc
      rw [hP, hP] at hc
      simpa using hc
    · rintro (⟨h0, h1⟩ | ⟨hcol, _⟩)
      · simp [h0, h1]
      · exact absurd hcol hD

/-- Shifting both arguments by lattice vectors leaves `Φ_q` of the difference unchanged. -/
theorem Φ_sub_congr (q : Model) {x x' y y' : ℝ} (hx : x ≡ x' [PMOD period q])
    (hy : y ≡ y' [PMOD period q]) : Φ q (x' - y') = Φ q (x - y) := by
  obtain ⟨z1, hz1⟩ := hx
  obtain ⟨z2, hz2⟩ := hy
  rw [show x' - y' = (x - y) + (z1 - z2) • period q by
    rw [sub_zsmul]; linarith]
  exact (Φ_periodic q).zsmul (z1 - z2) _

/-- **The exact fit has zero loss** (both models; directions modulo the kernel
period): `(θ_0,θ_1) ≡ (β_0,β_1)`, `c = s`, or swapped. -/
theorem L_eq_zero_of_exactFit (q : Model) (β s c θ : Fin 2 → ℝ)
    (h : (θ 0 ≡ β 0 [PMOD period q] ∧ θ 1 ≡ β 1 [PMOD period q] ∧ c 0 = s 0 ∧ c 1 = s 1) ∨
      (θ 0 ≡ β 1 [PMOD period q] ∧ θ 1 ≡ β 0 [PMOD period q] ∧ c 0 = s 1 ∧ c 1 = s 0)) :
    L q β s c θ = 0 := by
  unfold L kernelExcessLoss kernelMassLoss kernelTeacherSelfEnergy
  simp only [Fin.sum_univ_two]
  have hev : Φ q (β 1 - β 0) = Φ q (β 0 - β 1) := by rw [← Φ_even, neg_sub]
  have hr : ∀ x : ℝ, x ≡ x [PMOD period q] := fun x => AddCommGroup.modEq_refl x
  rcases h with ⟨h0, h1, hc0, hc1⟩ | ⟨h0, h1, hc0, hc1⟩
  · rw [← Φ_sub_congr q h0 h0, ← Φ_sub_congr q h0 h1, ← Φ_sub_congr q h1 h0,
      ← Φ_sub_congr q h1 h1, ← Φ_sub_congr q h0 (hr (β 0)), ← Φ_sub_congr q h0 (hr (β 1)),
      ← Φ_sub_congr q h1 (hr (β 0)), ← Φ_sub_congr q h1 (hr (β 1)), hc0, hc1]
    ring
  · rw [← Φ_sub_congr q h0 h0, ← Φ_sub_congr q h0 h1, ← Φ_sub_congr q h1 h0,
      ← Φ_sub_congr q h1 h1, ← Φ_sub_congr q h0 (hr (β 0)), ← Φ_sub_congr q h0 (hr (β 1)),
      ← Φ_sub_congr q h1 (hr (β 0)), ← Φ_sub_congr q h1 (hr (β 1)), hc0, hc1, hev]
    ring

/-- A zero student measure against a zero teacher field has zero loss. -/
theorem L_eq_zero_of_zeroMeasure (hP : ∀ t, P q β s t = 0)
    (hE : kernelTeacherSelfEnergy (Φ q) β s = 0)
    (hc : (c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD period q] ∧ c 0 + c 1 = 0)) :
    L q β s c θ = 0 := by
  have hcross : ∀ i, ∑ k, c i * s k * Φ q (θ i - β k) = c i * P q β s (θ i) := by
    intro i; unfold P; rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun k _ => by ring
  unfold L kernelExcessLoss kernelMassLoss
  rw [hE, Finset.sum_congr rfl fun i _ => hcross i]
  simp only [hP, mul_zero, Finset.sum_const_zero, sub_zero, add_zero, Fin.sum_univ_two]
  rcases hc with ⟨h0, h1⟩ | ⟨hcol, hsum⟩
  · simp [h0, h1]
  · have hev : Φ q (θ 1 - θ 0) = Φ q (θ 0 - θ 1) := by rw [← Φ_even, neg_sub]
    rw [Φ_of_modEq q hcol, hev, Φ_of_modEq q hcol, sub_self, Φ_zero, sub_self, Φ_zero,
      show c 1 = -c 0 by linarith]
    ring

end Rows

/-! ## 6. The collision block (`prop-collision-block`) -/

section Collision

variable {β s c θ : Fin 2 → ℝ} {q : Model}

/-- At collided students the loss depends on the masses only through their
sum, and reads the teacher potential at the common direction
(`eq-collision-total-mass-only`, hence items 1 and 2 of `prop-collision-block`:
flatness along the split line and the rigid-rotation formula). -/
theorem L_collided (q : Model) (β s c : Fin 2 → ℝ) (t : ℝ) :
    L q β s c ![t, t] =
      (1 / (2 * π)) * ((κ q / 2) * (c 0 + c 1) ^ 2 - (c 0 + c 1) * P q β s t +
        kernelTeacherSelfEnergy (Φ q) β s) := by
  unfold L kernelExcessLoss kernelMassLoss P
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    sub_self, Φ_zero]
  ring

/-- Item 1 of `prop-collision-block`: the split line is a level set. -/
theorem L_collided_split (q : Model) (β s : Fin 2 → ℝ) (a b a' b' t : ℝ) (h : a' + b' = a + b) :
    L q β s ![a', b'] ![t, t] = L q β s ![a, b] ![t, t] := by
  rw [L_collided, L_collided]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, h]

/-- The four rows at a collision: the mass pin `κ_q(c_0+c_1) = P_q(t)` and,
unless both students are dead, the torque root `A_q(t) = 0`. -/
theorem fourRows_collided_iff (q : Model) (β s c : Fin 2 → ℝ) (t : ℝ) :
    FourRows q β s c ![t, t] ↔
      κ q * (c 0 + c 1) = P q β s t ∧ ((c 0 ≠ 0 ∨ c 1 ≠ 0) → A q β s t = 0) := by
  unfold FourRows
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons, sub_self, Φ_zero,
    Hq_zero, mul_zero, sub_zero, add_zero]
  constructor
  · rintro ⟨_, h1, h2, h3⟩
    refine ⟨by linarith, fun hne => ?_⟩
    rcases hne with hne | hne
    · exact (mul_eq_zero.mp h2).resolve_left hne
    · exact (mul_eq_zero.mp h3).resolve_left hne
  · rintro ⟨hpin, hA⟩
    by_cases h0 : c 0 = 0
    · by_cases h1 : c 1 = 0
      · simp [h0, h1] at hpin ⊢
        linarith
      · have hAt := hA (Or.inr h1)
        exact ⟨by linarith, by linarith, by simp [h0], by simp [hAt]⟩
    · have hAt := hA (Or.inl h0)
      exact ⟨by linarith, by linarith, by simp [hAt], by simp [hAt]⟩

/-- The collision selector of the manuscript: the angular block `C_q` is
positive definite iff `τ_q P_q > 0 ∧ W c_0 c_1 τ_q < 0`
(`eq-collision-selector`, `clm:census-centered-selector`); on the degenerate
boundary `τ_q = 0` the flat selector `W'(t) = 0 ∧ c_0 c_1 < 0` decides, off
the kink lattice (`eq-flat-selector`, `lem-flat-sufficiency`). -/
def CollisionSelector (q : Model) (β s c : Fin 2 → ℝ) (t : ℝ) : Prop :=
  (0 < τ q β s t * P q β s t ∧ W β s t * (c 0 * c 1) * τ q β s t < 0) ∨
  (τ q β s t = 0 ∧ OffKink β t ∧ Wslope β s t = 0 ∧ c 0 * c 1 < 0)

/-- The angular block of the collision expansion (`eq-collision-angular-block`)
and its determinant/trace. -/
theorem collisionAngularBlock_det_trace (κq a b Pt Wt : ℝ) (hpin : κq * (a + b) = Pt) :
    (a * (Pt - 2 * Wt) - κq * a * b) * (b * (Pt - 2 * Wt) - κq * a * b) - (κq * a * b) ^ 2 =
        -2 * Wt * a * b * (Pt - 2 * Wt) ∧
      (a * (Pt - 2 * Wt) - κq * a * b) + (b * (Pt - 2 * Wt) - κq * a * b) =
        (Pt - 2 * Wt) * (a + b) - 2 * κq * a * b := by
  constructor
  · linear_combination (-(a * b * (Pt - 2 * Wt))) * hpin
  · ring

/-- `clm:census-collision-weight`: for a genuine teacher every torque root has
`W(t) ≠ 0`. -/
theorem W_ne_zero_of_torqueRoot_centered (h : CenteredGenuine β s) {t : ℝ}
    (hA : A .centered β s t = 0) : W β s t ≠ 0 := by
  sorry

theorem W_ne_zero_of_torqueRoot_plain (h : PlainGenuine β s) {t : ℝ}
    (hA : A .plainRelu β s t = 0) : W β s t ≠ 0 := by
  sorry

/-- `clm:census-centered-mixed-curvature`: for a mixed centered teacher every
torque root has `τ_C(t) ≠ 0`. -/
theorem τ_ne_zero_of_torqueRoot_mixed (h : CenteredGenuine β s) (hmix : s 0 * s 1 < 0)
    {t : ℝ} (hA : A .centered β s t = 0) : τ .centered β s t ≠ 0 := by
  sorry

/-- **The collision type** (`prop-collision-block`(4), `lem-flat-sufficiency`,
`stp:census-flat-*`, `stp:census-centered-collision`,
`stp:census-plain-gap-collision`): for a teacher with two genuine lines/rays, a
collided critical pair is a local minimum iff the collision selector holds. -/
theorem collided_isLocalMin_iff_selector (q : Model) (β s c : Fin 2 → ℝ) (t : ℝ)
    (hgen : (q = .centered ∧ CenteredGenuine β s) ∨ (q = .plainRelu ∧ PlainGenuine β s))
    (hA : A q β s t = 0) (hpin : κ q * (c 0 + c 1) = P q β s t) :
    IsLocalMinL q β s c ![t, t] ↔ CollisionSelector q β s c t := by
  sorry

/-- A collided pair against a teacher with two genuine lines/rays never has
zero loss (`prop-circle-l2`): the collision minima are spurious. -/
theorem collided_not_isGlobalMin (q : Model) (β s c : Fin 2 → ℝ) (t : ℝ)
    (hgen : (q = .centered ∧ CenteredGenuine β s) ∨ (q = .plainRelu ∧ PlainGenuine β s)) :
    ¬ IsGlobalMin q β s c ![t, t] := by
  sorry

end Collision

/-! ## 7. Dead students (`stp:census-centered-revival`, `clm:census-signed-revival`,
`clm:census-centered-three-translates`, `clm:census-plain-three-translates`) -/

section Dead

variable {β s c θ : Fin 2 → ℝ} {q : Model}

/-- `clm:census-signed-revival`: reviving a dead student at direction `x` with
mass `ε` changes the loss by exactly `(κ_q/2)ε² + ε F_q(x)` (times `1/(2π)`). -/
theorem L_revive_dead (q : Model) (β s c θ : Fin 2 → ℝ) (hdead : c 0 = 0) (ε x : ℝ) :
    L q β s ![ε, c 1] ![x, θ 1] - L q β s c θ =
      (1 / (2 * π)) * ((κ q / 2) * ε ^ 2 + ε * residual q β s c θ x) := by
  have h1 : Φ q (θ 1 - θ 1) = κ q := by rw [sub_self, Φ_zero]
  have hx : Φ q (x - x) = κ q := by rw [sub_self, Φ_zero]
  have hev : Φ q (θ 1 - x) = Φ q (x - θ 1) := by rw [← Φ_even, neg_sub]
  unfold L kernelExcessLoss kernelMassLoss residual P
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    h1, hx, hev, hdead]
  ring

/-- **Dead students are never local minima for a genuine teacher**
(`stp:census-centered-revival`, plain analogue in `stp:census-plain-gap-collision` /
`stp:census-plain-antipodal-separated`): the residual cannot vanish on an open
arc (three-translate rigidity), so a signed revival descends. -/
theorem dead_not_isLocalMin (q : Model) (β s c θ : Fin 2 → ℝ)
    (hgen : (q = .centered ∧ CenteredGenuine β s) ∨ (q = .plainRelu ∧ PlainTwoRay β s))
    (hdead : ∃ i, c i = 0) : ¬ IsLocalMinL q β s c θ := by
  sorry

/-- `clm:census-centered-three-translates`: no nontrivial combination of three
translates of `Φ_C` at distinct lines vanishes on a nonempty open arc. -/
theorem centered_three_translates_rigidity (a : Fin 3 → ℝ) (x : Fin 3 → ℝ)
    (hdist : ∀ i j, i ≠ j → ¬ x i ≡ x j [PMOD π]) (l r : ℝ) (hlr : l < r)
    (hzero : ∀ γ ∈ Ioo l r, ∑ i, a i * phiCos (γ - x i) = 0) : ∀ i, a i = 0 := by
  sorry

/-- `clm:census-plain-three-translates`. -/
theorem plain_three_translates_rigidity (a : Fin 3 → ℝ) (x : Fin 3 → ℝ)
    (hdist : ∀ i j, i ≠ j → ¬ x i ≡ x j [PMOD (2 * π)]) (l r : ℝ) (hlr : l < r)
    (hzero : ∀ γ ∈ Ioo l r, ∑ i, a i * phiCosJ (γ - x i) = 0) : ∀ i, a i = 0 := by
  sorry

end Dead

/-! ## 8. The centered families -/

section CenteredFamilies

/-- The exact fit `(θ_0,θ_1) ≡ (β_0,β_1)`, `(c_0,c_1) = (s_0,s_1)`, or swapped
(directions in `ℝ/πℤ`). -/
def CenteredExactFit (β s c θ : Fin 2 → ℝ) : Prop :=
  (θ 0 ≡ β 0 [PMOD π] ∧ θ 1 ≡ β 1 [PMOD π] ∧ c 0 = s 0 ∧ c 1 = s 1) ∨
  (θ 0 ≡ β 1 [PMOD π] ∧ θ 1 ≡ β 0 [PMOD π] ∧ c 0 = s 1 ∧ c 1 = s 0)

/-- Collision line: `θ_0 ≡ θ_1 ≡ t` a torque root, `κ_q(c_0+c_1) = P_q(t)`,
split free. -/
def CollisionLine (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  θ 0 ≡ θ 1 [PMOD period q] ∧ A q β s (θ 0) = 0 ∧ κ q * (c 0 + c 1) = P q β s (θ 0)

/-- One student dead, the other at a torque root with the solved mass
`κ_q c_j = P_q(θ_j) ≠ 0`; the dead direction solves the dead radial row
`c_j Φ_q(x − θ_j) = P_q(x)`. -/
def OneDead (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  ∃ i j : Fin 2, i ≠ j ∧ c i = 0 ∧ c j ≠ 0 ∧ A q β s (θ j) = 0 ∧
    κ q * c j = P q β s (θ j) ∧ c j * Φ q (θ i - θ j) = P q β s (θ i)

/-- Both students dead at zeros of the teacher potential. -/
def BothDead (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  c 0 = 0 ∧ c 1 = 0 ∧ P q β s (θ 0) = 0 ∧ P q β s (θ 1) = 0

/-- The bisector `(β_0+β_1)/2` of the two teacher directions. -/
def Bisector (β : Fin 2 → ℝ) : ℝ := (β 0 + β 1) / 2

/-- The bisector pair `{β̂/2, β̂/2 + π/2}` in `ℝ/πℤ`, in the original frame. -/
def QuarterPairAngles (β θ : Fin 2 → ℝ) : Prop :=
  (θ 0 ≡ Bisector β [PMOD π] ∧ θ 1 ≡ Bisector β + π / 2 [PMOD π]) ∨
  (θ 0 ≡ Bisector β + π / 2 [PMOD π] ∧ θ 1 ≡ Bisector β [PMOD π])

/-- Masses given by Cramér's rule (`eq-separated-cramer`). -/
def CramerMasses (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  c 0 = cramer0 q β s θ ∧ c 1 = cramer1 q β s θ

/-- Masses given by the angular rows (`eq-classification-angular-derivative-masses`),
in cleared form. -/
def TorqueMasses (q : Model) (β s c θ : Fin 2 → ℝ) : Prop :=
  c 1 * Hq q (θ 0 - θ 1) = A q β s (θ 0) ∧ c 0 * Hq q (θ 0 - θ 1) = - A q β s (θ 1)

/-- Strict interlacing of the two student lines with the two teacher lines in
`ℝ/πℤ`: in the frame `β_0 = 0`, one student in `(0, β̂)` and the other in
`(β̂, π)` (the beam coordinates `θ_1 = β̂ − a`, `θ_0 = β̂ + b`, or swapped). -/
def CenteredInterlaced (β θ : Fin 2 → ℝ) : Prop :=
  (0 < dirRep (θ 1 - β 0) ∧ dirRep (θ 1 - β 0) < canonicalGapC β ∧
      canonicalGapC β < dirRep (θ 0 - β 0)) ∨
  (0 < dirRep (θ 0 - β 0) ∧ dirRep (θ 0 - β 0) < canonicalGapC β ∧
      canonicalGapC β < dirRep (θ 1 - β 0))

/-- The orthogonal-diagonal common mass `√2 s_0(π+4)/(2(π+2))`. -/
def diagonalMass (s0 : ℝ) : ℝ := sqrt 2 * s0 * (π + 4) / (2 * (π + 2))

/-- The rows of the two centered tables. -/
inductive CenteredFamily where
  | zeroMeasure
  | oneLineExact
  | oneLinePerp
  | exactFit
  | collision
  | quarterPair
  | orthogonalDiagonal
  | separatedBeam
  | oneDead
  | bothDead
  deriving DecidableEq

/-- Membership in each row of the centered tables (`thm-centered-two-student-enumeration`). -/
def CenteredMember : CenteredFamily → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → Prop
  | .zeroMeasure, β, s, c, θ =>
      CenteredZeroField β s ∧
        ((c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD π] ∧ c 0 + c 1 = 0))
  | .oneLineExact, β, s, c, θ =>
      ∃ a σ, CenteredOneLine β s a σ ∧ (∀ i, c i = 0 ∨ θ i ≡ a [PMOD π]) ∧ c 0 + c 1 = σ
  | .oneLinePerp, β, s, c, θ =>
      ∃ a σ, CenteredOneLine β s a σ ∧ θ 0 ≡ a + π / 2 [PMOD π] ∧ θ 1 ≡ a + π / 2 [PMOD π] ∧
        c 0 + c 1 = 2 * σ / π
  | .exactFit, β, s, c, θ => CenteredGenuine β s ∧ CenteredExactFit β s c θ
  | .collision, β, s, c, θ => CenteredGenuine β s ∧ CollisionLine .centered β s c θ
  | .quarterPair, β, s, c, θ =>
      CenteredGenuine β s ∧ s 0 = s 1 ∧ ¬ β 1 - β 0 ≡ π / 2 [PMOD π] ∧
        QuarterPairAngles β θ ∧ CramerMasses .centered β s c θ
  | .orthogonalDiagonal, β, s, c, θ =>
      CenteredGenuine β s ∧ s 0 = s 1 ∧ β 1 - β 0 ≡ π / 2 [PMOD π] ∧
        QuarterPairAngles β θ ∧ c 0 = diagonalMass (s 0) ∧ c 1 = diagonalMass (s 0)
  | .separatedBeam, β, s, c, θ =>
      CenteredGenuine β s ∧ 0 < s 0 * s 1 ∧ s 0 ≠ s 1 ∧ c 0 ≠ 0 ∧ c 1 ≠ 0 ∧
        CenteredInterlaced β θ ∧ CramerMasses .centered β s c θ ∧ TorqueMasses .centered β s c θ
  | .oneDead, β, s, c, θ => CenteredGenuine β s ∧ OneDead .centered β s c θ
  | .bothDead, β, s, c, θ => CenteredGenuine β s ∧ BothDead .centered β s c θ

/-- The union of all rows. -/
def CenteredFamilies (β s c θ : Fin 2 → ℝ) : Prop := ∃ f, CenteredMember f β s c θ

/-- The variational label each row carries. -/
def centeredFamilyType : CenteredFamily → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → Prop
  | .zeroMeasure, β, s, c, θ => IsGlobalMin .centered β s c θ
  | .oneLineExact, β, s, c, θ => IsGlobalMin .centered β s c θ
  | .oneLinePerp, β, s, c, θ => IsTopologicalSaddle .centered β s c θ
  | .exactFit, β, s, c, θ => IsGlobalMin .centered β s c θ
  | .collision, β, s, c, θ =>
      (CollisionSelector .centered β s c (θ 0) ∧ IsSpuriousMin .centered β s c θ) ∨
      (¬ CollisionSelector .centered β s c (θ 0) ∧ IsTopologicalSaddle .centered β s c θ)
  | .quarterPair, β, s, c, θ => IsTopologicalSaddle .centered β s c θ
  | .orthogonalDiagonal, β, s, c, θ => IsTopologicalSaddle .centered β s c θ
  | .separatedBeam, β, s, c, θ => IsTopologicalSaddle .centered β s c θ
  | .oneDead, β, s, c, θ => IsTopologicalSaddle .centered β s c θ
  | .bothDead, β, s, c, θ => IsTopologicalSaddle .centered β s c θ

/-! ### Per-stratum classification lemmas -/

/-- `stp:census-centered-one-line` (zero-field rows): for a zero teacher field,
critical points are exactly the zero student measures. -/
theorem centered_zeroField_critical_iff {β s c θ : Fin 2 → ℝ} (h : CenteredZeroField β s) :
    IsCritical .centered β s c θ ↔
      ((c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD π] ∧ c 0 + c 1 = 0)) :=
  zeroField_critical_iff (q := .centered) (P_eq_zero_of_centeredZeroField h)
    (A_eq_zero_of_centeredZeroField h)

/-- `stp:census-centered-one-line` (one-line rows): the exact line (global) and
the perpendicular line (saddle); no separated two-mass point exists. -/
theorem centered_oneLine_critical_iff {β s c θ : Fin 2 → ℝ} {a σ : ℝ}
    (h : CenteredOneLine β s a σ) :
    IsCritical .centered β s c θ ↔
      (((∀ i, c i = 0 ∨ θ i ≡ a [PMOD π]) ∧ c 0 + c 1 = σ) ∨
        (θ 0 ≡ a + π / 2 [PMOD π] ∧ θ 1 ≡ a + π / 2 [PMOD π] ∧ c 0 + c 1 = 2 * σ / π)) := by
  sorry

/-- `stp:census-centered-zero-mass`, `stp:census-centered-quarter`,
`stp:census-centered-beam`, `stp:census-centered-mixed`: the genuine stratum. -/
theorem centered_genuine_critical_iff {β s c θ : Fin 2 → ℝ} (h : CenteredGenuine β s) :
    IsCritical .centered β s c θ ↔
      (CenteredExactFit β s c θ ∨ CollisionLine .centered β s c θ ∨
        (s 0 = s 1 ∧ ¬ β 1 - β 0 ≡ π / 2 [PMOD π] ∧ QuarterPairAngles β θ ∧
          CramerMasses .centered β s c θ) ∨
        (s 0 = s 1 ∧ β 1 - β 0 ≡ π / 2 [PMOD π] ∧ QuarterPairAngles β θ ∧
          c 0 = diagonalMass (s 0) ∧ c 1 = diagonalMass (s 0)) ∨
        (0 < s 0 * s 1 ∧ s 0 ≠ s 1 ∧ c 0 ≠ 0 ∧ c 1 ≠ 0 ∧ CenteredInterlaced β θ ∧
          CramerMasses .centered β s c θ ∧ TorqueMasses .centered β s c θ) ∨
        OneDead .centered β s c θ ∨ BothDead .centered β s c θ) := by
  sorry

/-- **Complete centered two-student enumeration** (`thm:centered-two-teacher-census`,
long tree `centered_two_student_full_classification`): a tuple is critical for
`L_C` iff it belongs to one of the ten rows.  Assembled from the three
stratum lemmas by the exhaustive trichotomy. -/
theorem centered_two_student_classification (β s c θ : Fin 2 → ℝ) :
    IsCriticalL_C β s c θ ↔ CenteredFamilies β s c θ := by
  unfold CenteredFamilies
  constructor
  · intro hcrit
    rcases centered_strata_exhaustive β s with hz | ⟨a, σ, hone⟩ | hgen
    · exact ⟨.zeroMeasure, hz, (centered_zeroField_critical_iff hz).mp hcrit⟩
    · rcases (centered_oneLine_critical_iff hone).mp hcrit with hex | hperp
      · exact ⟨.oneLineExact, a, σ, hone, hex.1, hex.2⟩
      · exact ⟨.oneLinePerp, a, σ, hone, hperp.1, hperp.2.1, hperp.2.2⟩
    · rcases (centered_genuine_critical_iff hgen).mp hcrit with
        h | h | h | h | h | h | h
      · exact ⟨.exactFit, hgen, h⟩
      · exact ⟨.collision, hgen, h⟩
      · exact ⟨.quarterPair, hgen, h⟩
      · exact ⟨.orthogonalDiagonal, hgen, h⟩
      · exact ⟨.separatedBeam, hgen, h⟩
      · exact ⟨.oneDead, hgen, h⟩
      · exact ⟨.bothDead, hgen, h⟩
  · rintro ⟨f, hf⟩
    cases f with
    | zeroMeasure => exact (centered_zeroField_critical_iff hf.1).mpr hf.2
    | oneLineExact =>
        obtain ⟨a, σ, hone, hc, hsum⟩ := hf
        exact (centered_oneLine_critical_iff hone).mpr (Or.inl ⟨hc, hsum⟩)
    | oneLinePerp =>
        obtain ⟨a, σ, hone, h0, h1, hsum⟩ := hf
        exact (centered_oneLine_critical_iff hone).mpr (Or.inr ⟨h0, h1, hsum⟩)
    | exactFit => exact (centered_genuine_critical_iff hf.1).mpr (Or.inl hf.2)
    | collision => exact (centered_genuine_critical_iff hf.1).mpr (Or.inr (Or.inl hf.2))
    | quarterPair =>
        exact (centered_genuine_critical_iff hf.1).mpr (Or.inr (Or.inr (Or.inl hf.2)))
    | orthogonalDiagonal =>
        exact (centered_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))
    | separatedBeam =>
        exact (centered_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2)))))
    | oneDead =>
        exact (centered_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))))
    | bothDead =>
        exact (centered_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr hf.2))))))

/-! ### Per-family type lemmas -/

/-- The zero student measure against a zero field has zero loss. -/
theorem centered_zeroMeasure_global {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .zeroMeasure β s c θ) : IsGlobalMin .centered β s c θ :=
  isGlobalMin_of_L_eq_zero (L_eq_zero_of_zeroMeasure (P_eq_zero_of_centeredZeroField h.1)
    (selfEnergy_eq_zero_of_centeredZeroField h.1) h.2)

/-- The exact one-line fit has zero loss. -/
theorem centered_oneLineExact_global {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .oneLineExact β s c θ) : IsGlobalMin .centered β s c θ := by
  sorry

/-- The perpendicular line is a saddle: co-rotating the pair decreases the loss
at second order, `−(2/π)σ²` (`stp:census-centered-one-line`). -/
theorem centered_oneLinePerp_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .oneLinePerp β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  sorry

/-- The exact fit has zero loss. -/
theorem centered_exactFit_global {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .exactFit β s c θ) : IsGlobalMin .centered β s c θ :=
  isGlobalMin_of_L_eq_zero (L_eq_zero_of_exactFit .centered β s c θ h.2)

/-- Collision rows: selector decides (`stp:census-centered-collision`). -/
theorem centered_collision_type {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .collision β s c θ) : centeredFamilyType .collision β s c θ := by
  sorry

/-- `stp:census-centered-quarter`: the quarter pair is a saddle. -/
theorem centered_quarterPair_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .quarterPair β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  sorry

/-- `stp:census-centered-quarter`: the orthogonal diagonal is a saddle. -/
theorem centered_orthogonalDiagonal_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .orthogonalDiagonal β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  sorry

/-- `stp:census-centered-beam`, `clm:census-beam-signs`: the separated beam
root is a saddle (the hard separated second-variation step). -/
theorem centered_separatedBeam_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .separatedBeam β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  sorry

/-- `stp:census-centered-revival`: a dead student against a genuine teacher is a saddle. -/
theorem centered_oneDead_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .oneDead β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  have hcrit : IsCritical .centered β s c θ :=
    (centered_two_student_classification β s c θ).mpr ⟨.oneDead, h⟩
  rw [isTopologicalSaddle_iff_not_isLocalMin hcrit]
  obtain ⟨i, j, _, hi, _⟩ := h.2
  exact dead_not_isLocalMin .centered β s c θ (Or.inl ⟨rfl, h.1⟩) ⟨i, hi⟩

theorem centered_bothDead_saddle {β s c θ : Fin 2 → ℝ}
    (h : CenteredMember .bothDead β s c θ) : IsTopologicalSaddle .centered β s c θ := by
  have hcrit : IsCritical .centered β s c θ :=
    (centered_two_student_classification β s c θ).mpr ⟨.bothDead, h⟩
  rw [isTopologicalSaddle_iff_not_isLocalMin hcrit]
  exact dead_not_isLocalMin .centered β s c θ (Or.inl ⟨rfl, h.1⟩) ⟨0, h.2.1⟩

/-- **Each centered row carries its displayed label**
(`thm:centered-two-teacher-census`, type clause; long tree
`centered_two_student_variational_type`). -/
theorem centered_two_student_type (f : CenteredFamily) (β s c θ : Fin 2 → ℝ)
    (h : CenteredMember f β s c θ) : centeredFamilyType f β s c θ := by
  cases f with
  | zeroMeasure => exact centered_zeroMeasure_global h
  | oneLineExact => exact centered_oneLineExact_global h
  | oneLinePerp => exact centered_oneLinePerp_saddle h
  | exactFit => exact centered_exactFit_global h
  | collision => exact centered_collision_type h
  | quarterPair => exact centered_quarterPair_saddle h
  | orthogonalDiagonal => exact centered_orthogonalDiagonal_saddle h
  | separatedBeam => exact centered_separatedBeam_saddle h
  | oneDead => exact centered_oneDead_saddle h
  | bothDead => exact centered_bothDead_saddle h

/-- Every centered critical point is a global minimum, a spurious local
minimum, or a topological saddle. -/
theorem centered_variationalType {β s c θ : Fin 2 → ℝ} (h : IsCriticalL_C β s c θ) :
    VariationalType .centered β s c θ :=
  variationalType_of_isCritical h

/-! ### Structure of the separated beam stratum -/

/-- `stp:census-centered-beam` (structure clause): for a same-sign teacher with
`s_0 ≠ s_1`, the beam row has positive student masses and exists at exactly
one interlaced angle pair up to the student swap
(`centeredBeamRoot_existsUnique` in the long tree). -/
theorem centered_separatedBeam_existsUnique {β s : Fin 2 → ℝ} (hgen : CenteredGenuine β s)
    (hsame : 0 < s 0 * s 1) (hne : s 0 ≠ s 1) :
    (∃ c θ : Fin 2 → ℝ, CenteredMember .separatedBeam β s c θ) ∧
    (∀ c θ c' θ' : Fin 2 → ℝ, CenteredMember .separatedBeam β s c θ →
      CenteredMember .separatedBeam β s c' θ' →
      ((θ 0 ≡ θ' 0 [PMOD π] ∧ θ 1 ≡ θ' 1 [PMOD π]) ∨
        (θ 0 ≡ θ' 1 [PMOD π] ∧ θ 1 ≡ θ' 0 [PMOD π]))) ∧
    (∀ c θ : Fin 2 → ℝ, CenteredMember .separatedBeam β s c θ → 0 < c 0 * s 0 ∧ 0 < c 1 * s 0) := by
  sorry

/-- `clm:census-beam-signs`: `G(a,b) > 0` and `Θ(a,b) > 0` on the open beam
triangle `a, b > 0`, `a + b < π`. -/
theorem beam_signs {a b : ℝ} (ha : 0 < a) (hb : 0 < b) (hab : a + b < π) :
    0 < (a + b) * b * sin a - a * sin (a + b) * sin b ∧
    0 < ((a + b) ^ 2 - sin (a + b) ^ 2) * (sin a - a * cos a) * (b * sin b) -
      (a * sin b - b * sin a) * (a * sin b + b * sin a) *
        (sin (a + b) - (a + b) * cos (a + b)) := by
  sorry

end CenteredFamilies

/-! ## 9. The plain-ReLU families -/

section PlainFamilies

/-- The plain exact fit, directions in `ℝ/2πℤ`.  (For a two-ray teacher this is
the manuscript's "`(0,β̂)` in `ℝ/πℤ` with `u = v`".) -/
def PlainExactFit (β s c θ : Fin 2 → ℝ) : Prop :=
  (θ 0 ≡ β 0 [PMOD (2 * π)] ∧ θ 1 ≡ β 1 [PMOD (2 * π)] ∧ c 0 = s 0 ∧ c 1 = s 1) ∨
  (θ 0 ≡ β 1 [PMOD (2 * π)] ∧ θ 1 ≡ β 0 [PMOD (2 * π)] ∧ c 0 = s 1 ∧ c 1 = s 0)

/-- The endpoint null cubic `N(a,b) = ab|μ|³ − S(a|b|³ + b|a|³)`
(`eq-antipodal-endpoint-cusp-selector`), with `μ = a + b`. -/
def nullCubic (S a b : ℝ) : ℝ := a * b * |a + b| ^ 3 - S * (a * |b| ^ 3 + b * |a| ^ 3)

/-- The raw angular block `𝒜` and the cleared angular Schur complement `T`
(`eq-separated-schur`) of a separated plain-ReLU pair. -/
def schurA00 (β s c θ : Fin 2 → ℝ) : ℝ :=
  c 0 * τ .plainRelu β s (θ 0) - c 0 * c 1 * Hslope .plainRelu (θ 0 - θ 1)
def schurA11 (β s c θ : Fin 2 → ℝ) : ℝ :=
  c 1 * τ .plainRelu β s (θ 1) - c 0 * c 1 * Hslope .plainRelu (θ 0 - θ 1)
def schurA01 (c θ : Fin 2 → ℝ) : ℝ := c 0 * c 1 * Hslope .plainRelu (θ 0 - θ 1)
def schurT00 (β s c θ : Fin 2 → ℝ) : ℝ :=
  Δ .plainRelu (θ 0 - θ 1) * schurA00 β s c θ - π * Hq .plainRelu (θ 0 - θ 1) ^ 2 * c 0 ^ 2
def schurT11 (β s c θ : Fin 2 → ℝ) : ℝ :=
  Δ .plainRelu (θ 0 - θ 1) * schurA11 β s c θ - π * Hq .plainRelu (θ 0 - θ 1) ^ 2 * c 1 ^ 2
def schurT01 (c θ : Fin 2 → ℝ) : ℝ :=
  Δ .plainRelu (θ 0 - θ 1) * schurA01 c θ -
    Hq .plainRelu (θ 0 - θ 1) ^ 2 * Φ .plainRelu (θ 0 - θ 1) * c 0 * c 1
def schurDet (β s c θ : Fin 2 → ℝ) : ℝ :=
  schurT00 β s c θ * schurT11 β s c θ - schurT01 c θ ^ 2

/-- The Schur-typed label of a smooth separated plain-ReLU pair. -/
def SchurType (β s c θ : Fin 2 → ℝ) : Prop :=
  (0 < schurT00 β s c θ ∧ 0 < schurDet β s c θ → IsSpuriousMin .plainRelu β s c θ) ∧
  (schurDet β s c θ < 0 ∨ schurT00 β s c θ < 0 ∨ schurT11 β s c θ < 0 →
    IsTopologicalSaddle .plainRelu β s c θ) ∧
  (0 ≤ schurT00 β s c θ ∧ 0 ≤ schurT11 β s c θ ∧ schurDet β s c θ = 0 →
    IsSpuriousMin .plainRelu β s c θ ∨ IsTopologicalSaddle .plainRelu β s c θ)

/-- A student sits on a ray opposite a teacher ray (a ghost corner). -/
def GhostCorner (β θ : Fin 2 → ℝ) : Prop := ∃ i k, θ i ≡ β k + π [PMOD (2 * π)]

/-- No student on a teacher ray or on an opposite ray. -/
def SmoothPair (β θ : Fin 2 → ℝ) : Prop :=
  ∀ i k, ¬ θ i ≡ β k [PMOD (2 * π)] ∧ ¬ θ i ≡ β k + π [PMOD (2 * π)]

/-- Two live students at a gap `D ≢ 0, π`. -/
def PlainSeparatedLive (c θ : Fin 2 → ℝ) : Prop :=
  c 0 ≠ 0 ∧ c 1 ≠ 0 ∧ ¬ θ 0 ≡ θ 1 [PMOD (2 * π)] ∧ ¬ θ 0 ≡ θ 1 + π [PMOD (2 * π)]

/-- The symmetric separated pair of the antipodal stratum
(`eq-antipodal-midpoint-root`): common mass `a`, opening `E`, axis `φ`. -/
def AntipodalMidpointPair (β s c θ : Fin 2 → ℝ) : Prop :=
  ∃ φ E a : ℝ, 0 < E ∧ E < π ∧
    a = (s 0 + s 1) * sin (E / 2) / (E + sin E) ∧
    ((φ ≡ β 0 [PMOD (2 * π)] ∧
        π * s 0 = (s 0 + s 1) * (E / 2) + 2 * a * (π - E) * cos (E / 2)) ∨
      (φ ≡ β 1 [PMOD (2 * π)] ∧
        π * s 1 = (s 0 + s 1) * (E / 2) + 2 * a * (π - E) * cos (E / 2))) ∧
    c 0 = a ∧ c 1 = a ∧
    ((θ 0 ≡ φ + E / 2 [PMOD (2 * π)] ∧ θ 1 ≡ φ - E / 2 [PMOD (2 * π)]) ∨
      (θ 0 ≡ φ - E / 2 [PMOD (2 * π)] ∧ θ 1 ≡ φ + E / 2 [PMOD (2 * π)]))

/-- The rows of the three plain-ReLU tables. -/
inductive PlainFamily where
  | zeroMeasure
  | oneRayExact
  | oneRayAntipodalCancel
  | antipodalInteriorCollision
  | antipodalEndpointCollision
  | oppositePerp
  | separatedMidpointPair
  | exactFit
  | collisionRay
  | oppositeBisector
  | oneDead
  | bothDead
  | separatedGhost
  | separatedSmooth
  deriving DecidableEq

/-- Membership in each row of the plain-ReLU tables (`thm-plain-critical-enumeration`). -/
def PlainMember : PlainFamily → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → Prop
  | .zeroMeasure, β, s, c, θ =>
      PlainZeroField β s ∧
        ((c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD (2 * π)] ∧ c 0 + c 1 = 0))
  | .oneRayExact, β, s, c, θ =>
      ∃ a σ, PlainOneRay β s a σ ∧ (∀ i, c i = 0 ∨ θ i ≡ a [PMOD (2 * π)]) ∧ c 0 + c 1 = σ
  | .oneRayAntipodalCancel, β, s, c, θ =>
      ∃ a σ, PlainOneRay β s a σ ∧ θ 0 ≡ a + π [PMOD (2 * π)] ∧ θ 1 ≡ a + π [PMOD (2 * π)] ∧
        c 0 + c 1 = 0
  | .antipodalInteriorCollision, β, s, c, θ =>
      PlainAntipodal β s ∧ θ 0 ≡ θ 1 [PMOD (2 * π)] ∧
        ¬ θ 0 ≡ β 0 [PMOD (2 * π)] ∧ ¬ θ 0 ≡ β 1 [PMOD (2 * π)] ∧
        A .plainRelu β s (θ 0) = 0 ∧ π * (c 0 + c 1) = P .plainRelu β s (θ 0)
  | .antipodalEndpointCollision, β, s, c, θ =>
      PlainAntipodal β s ∧ ∃ k : Fin 2, θ 0 ≡ β k [PMOD (2 * π)] ∧ θ 1 ≡ β k [PMOD (2 * π)] ∧
        c 0 + c 1 = s k
  | .oppositePerp, β, s, c, θ =>
      PlainAntipodal β s ∧ s 0 = s 1 ∧
        ((θ 0 ≡ β 0 + π / 2 [PMOD (2 * π)] ∧ θ 1 ≡ β 0 + 3 * π / 2 [PMOD (2 * π)]) ∨
          (θ 0 ≡ β 0 + 3 * π / 2 [PMOD (2 * π)] ∧ θ 1 ≡ β 0 + π / 2 [PMOD (2 * π)])) ∧
        c 0 = 2 * s 0 / π ∧ c 1 = 2 * s 0 / π
  | .separatedMidpointPair, β, s, c, θ => PlainAntipodal β s ∧ AntipodalMidpointPair β s c θ
  | .exactFit, β, s, c, θ => PlainTwoRay β s ∧ PlainExactFit β s c θ
  | .collisionRay, β, s, c, θ => PlainGenuine β s ∧ CollisionLine .plainRelu β s c θ
  | .oppositeBisector, β, s, c, θ =>
      PlainGenuine β s ∧ s 0 = s 1 ∧ θ 0 ≡ Bisector β [PMOD π] ∧
        θ 1 ≡ θ 0 + π [PMOD (2 * π)] ∧
        c 0 = P .plainRelu β s (θ 0) / π ∧ c 1 = P .plainRelu β s (θ 1) / π
  | .oneDead, β, s, c, θ => PlainTwoRay β s ∧ OneDead .plainRelu β s c θ
  | .bothDead, β, s, c, θ => PlainTwoRay β s ∧ BothDead .plainRelu β s c θ
  | .separatedGhost, β, s, c, θ =>
      PlainGenuine β s ∧ PlainSeparatedLive c θ ∧ GhostCorner β θ ∧
        CramerMasses .plainRelu β s c θ ∧ TorqueMasses .plainRelu β s c θ
  | .separatedSmooth, β, s, c, θ =>
      PlainGenuine β s ∧ PlainSeparatedLive c θ ∧ SmoothPair β θ ∧
        CramerMasses .plainRelu β s c θ ∧ TorqueMasses .plainRelu β s c θ

/-- The union of all plain-ReLU rows. -/
def PlainFamilies (β s c θ : Fin 2 → ℝ) : Prop := ∃ f, PlainMember f β s c θ

/-- The variational label each plain-ReLU row carries. -/
def plainFamilyType : PlainFamily → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → (Fin 2 → ℝ) → Prop
  | .zeroMeasure, β, s, c, θ => IsGlobalMin .plainRelu β s c θ
  | .oneRayExact, β, s, c, θ => IsGlobalMin .plainRelu β s c θ
  | .oneRayAntipodalCancel, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .antipodalInteriorCollision, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .antipodalEndpointCollision, β, s, c, θ =>
      (0 < nullCubic (s 0 + s 1) (c 0) (c 1) ∧ IsSpuriousMin .plainRelu β s c θ) ∨
      (nullCubic (s 0 + s 1) (c 0) (c 1) ≤ 0 ∧ IsTopologicalSaddle .plainRelu β s c θ)
  | .oppositePerp, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .separatedMidpointPair, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .exactFit, β, s, c, θ => IsGlobalMin .plainRelu β s c θ
  | .collisionRay, β, s, c, θ =>
      (CollisionSelector .plainRelu β s c (θ 0) ∧ IsSpuriousMin .plainRelu β s c θ) ∨
      (¬ CollisionSelector .plainRelu β s c (θ 0) ∧ IsTopologicalSaddle .plainRelu β s c θ)
  | .oppositeBisector, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .oneDead, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .bothDead, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .separatedGhost, β, s, c, θ => IsTopologicalSaddle .plainRelu β s c θ
  | .separatedSmooth, β, s, c, θ => SchurType β s c θ

/-! ### Per-stratum classification lemmas -/

/-- `stp:census-plain-reduction` (zero-field rows). -/
theorem plain_zeroField_critical_iff {β s c θ : Fin 2 → ℝ} (h : PlainZeroField β s) :
    IsCritical .plainRelu β s c θ ↔
      ((c 0 = 0 ∧ c 1 = 0) ∨ (θ 0 ≡ θ 1 [PMOD (2 * π)] ∧ c 0 + c 1 = 0)) :=
  zeroField_critical_iff (q := .plainRelu) (P_eq_zero_of_plainZeroField h)
    (A_eq_zero_of_plainZeroField h)

/-- `stp:census-plain-reduction` (one-ray rows): the exact ray and the
antipodal cancellation. -/
theorem plain_oneRay_critical_iff {β s c θ : Fin 2 → ℝ} {a σ : ℝ} (h : PlainOneRay β s a σ) :
    IsCritical .plainRelu β s c θ ↔
      (((∀ i, c i = 0 ∨ θ i ≡ a [PMOD (2 * π)]) ∧ c 0 + c 1 = σ) ∨
        (θ 0 ≡ a + π [PMOD (2 * π)] ∧ θ 1 ≡ a + π [PMOD (2 * π)] ∧ c 0 + c 1 = 0)) := by
  sorry

/-- `stp:census-plain-antipodal-collision`, `stp:census-plain-antipodal-separated`:
the antipodal stratum. -/
theorem plain_antipodal_critical_iff {β s c θ : Fin 2 → ℝ} (h : PlainAntipodal β s) :
    IsCritical .plainRelu β s c θ ↔
      ((θ 0 ≡ θ 1 [PMOD (2 * π)] ∧ ¬ θ 0 ≡ β 0 [PMOD (2 * π)] ∧ ¬ θ 0 ≡ β 1 [PMOD (2 * π)] ∧
          A .plainRelu β s (θ 0) = 0 ∧ π * (c 0 + c 1) = P .plainRelu β s (θ 0)) ∨
        (∃ k : Fin 2, θ 0 ≡ β k [PMOD (2 * π)] ∧ θ 1 ≡ β k [PMOD (2 * π)] ∧ c 0 + c 1 = s k) ∨
        (s 0 = s 1 ∧
          ((θ 0 ≡ β 0 + π / 2 [PMOD (2 * π)] ∧ θ 1 ≡ β 0 + 3 * π / 2 [PMOD (2 * π)]) ∨
            (θ 0 ≡ β 0 + 3 * π / 2 [PMOD (2 * π)] ∧ θ 1 ≡ β 0 + π / 2 [PMOD (2 * π)])) ∧
          c 0 = 2 * s 0 / π ∧ c 1 = 2 * s 0 / π) ∨
        AntipodalMidpointPair β s c θ ∨
        PlainExactFit β s c θ ∨
        OneDead .plainRelu β s c θ ∨ BothDead .plainRelu β s c θ) := by
  sorry

/-- `stp:census-plain-gap-collision`, `stp:census-plain-gap-separated`,
`clm:census-plain-load-sign`, `clm:census-plain-ghost-diagonal`: the genuine
stratum.  The "real corner" row of the table is absorbed: a separated live
pair with a student on a teacher ray is the exact fit. -/
theorem plain_genuine_critical_iff {β s c θ : Fin 2 → ℝ} (h : PlainGenuine β s) :
    IsCritical .plainRelu β s c θ ↔
      (CollisionLine .plainRelu β s c θ ∨
        (s 0 = s 1 ∧ θ 0 ≡ Bisector β [PMOD π] ∧ θ 1 ≡ θ 0 + π [PMOD (2 * π)] ∧
          c 0 = P .plainRelu β s (θ 0) / π ∧ c 1 = P .plainRelu β s (θ 1) / π) ∨
        PlainExactFit β s c θ ∨
        OneDead .plainRelu β s c θ ∨ BothDead .plainRelu β s c θ ∨
        (PlainSeparatedLive c θ ∧ GhostCorner β θ ∧
          CramerMasses .plainRelu β s c θ ∧ TorqueMasses .plainRelu β s c θ) ∨
        (PlainSeparatedLive c θ ∧ SmoothPair β θ ∧
          CramerMasses .plainRelu β s c θ ∧ TorqueMasses .plainRelu β s c θ)) := by
  sorry

/-- The real-corner row (`stp:census-plain-gap-separated`): a separated live
critical pair with a student on a teacher ray is the exact fit. -/
theorem plain_realCorner_forces_exactFit {β s c θ : Fin 2 → ℝ} (h : PlainGenuine β s)
    (hcrit : IsCritical .plainRelu β s c θ) (hlive : PlainSeparatedLive c θ)
    (hcorner : ∃ i k, θ i ≡ β k [PMOD (2 * π)]) : PlainExactFit β s c θ := by
  sorry

/-- **Complete plain-ReLU two-student enumeration** (`thm:plain-two-teacher-census`,
long tree `noncentered_two_student_full_classification`).  Assembled from the
four stratum lemmas. -/
theorem plain_two_student_classification (β s c θ : Fin 2 → ℝ) :
    IsCriticalL_R β s c θ ↔ PlainFamilies β s c θ := by
  unfold PlainFamilies
  constructor
  · intro hcrit
    rcases plain_strata_exhaustive β s with hz | ⟨a, σ, hone⟩ | hanti | hgen
    · exact ⟨.zeroMeasure, hz, (plain_zeroField_critical_iff hz).mp hcrit⟩
    · rcases (plain_oneRay_critical_iff hone).mp hcrit with hex | hcancel
      · exact ⟨.oneRayExact, a, σ, hone, hex.1, hex.2⟩
      · exact ⟨.oneRayAntipodalCancel, a, σ, hone, hcancel.1, hcancel.2.1, hcancel.2.2⟩
    · rcases (plain_antipodal_critical_iff hanti).mp hcrit with
        h | h | h | h | h | h | h
      · exact ⟨.antipodalInteriorCollision, hanti, h⟩
      · exact ⟨.antipodalEndpointCollision, hanti, h⟩
      · exact ⟨.oppositePerp, hanti, h⟩
      · exact ⟨.separatedMidpointPair, hanti, h⟩
      · exact ⟨.exactFit, Or.inl hanti, h⟩
      · exact ⟨.oneDead, Or.inl hanti, h⟩
      · exact ⟨.bothDead, Or.inl hanti, h⟩
    · rcases (plain_genuine_critical_iff hgen).mp hcrit with
        h | h | h | h | h | h | h
      · exact ⟨.collisionRay, hgen, h⟩
      · exact ⟨.oppositeBisector, hgen, h⟩
      · exact ⟨.exactFit, Or.inr hgen, h⟩
      · exact ⟨.oneDead, Or.inr hgen, h⟩
      · exact ⟨.bothDead, Or.inr hgen, h⟩
      · exact ⟨.separatedGhost, hgen, h⟩
      · exact ⟨.separatedSmooth, hgen, h⟩
  · rintro ⟨f, hf⟩
    cases f with
    | zeroMeasure => exact (plain_zeroField_critical_iff hf.1).mpr hf.2
    | oneRayExact =>
        obtain ⟨a, σ, hone, hc, hsum⟩ := hf
        exact (plain_oneRay_critical_iff hone).mpr (Or.inl ⟨hc, hsum⟩)
    | oneRayAntipodalCancel =>
        obtain ⟨a, σ, hone, h0, h1, hsum⟩ := hf
        exact (plain_oneRay_critical_iff hone).mpr (Or.inr ⟨h0, h1, hsum⟩)
    | antipodalInteriorCollision =>
        exact (plain_antipodal_critical_iff hf.1).mpr (Or.inl hf.2)
    | antipodalEndpointCollision =>
        exact (plain_antipodal_critical_iff hf.1).mpr (Or.inr (Or.inl hf.2))
    | oppositePerp =>
        exact (plain_antipodal_critical_iff hf.1).mpr (Or.inr (Or.inr (Or.inl hf.2)))
    | separatedMidpointPair =>
        exact (plain_antipodal_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))
    | exactFit =>
        rcases hf.1 with hanti | hgen
        · exact (plain_antipodal_critical_iff hanti).mpr
            (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2)))))
        · exact (plain_genuine_critical_iff hgen).mpr (Or.inr (Or.inr (Or.inl hf.2)))
    | collisionRay => exact (plain_genuine_critical_iff hf.1).mpr (Or.inl hf.2)
    | oppositeBisector => exact (plain_genuine_critical_iff hf.1).mpr (Or.inr (Or.inl hf.2))
    | oneDead =>
        rcases hf.1 with hanti | hgen
        · exact (plain_antipodal_critical_iff hanti).mpr
            (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))))
        · exact (plain_genuine_critical_iff hgen).mpr
            (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))
    | bothDead =>
        rcases hf.1 with hanti | hgen
        · exact (plain_antipodal_critical_iff hanti).mpr
            (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr hf.2))))))
        · exact (plain_genuine_critical_iff hgen).mpr
            (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2)))))
    | separatedGhost =>
        exact (plain_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inl hf.2))))))
    | separatedSmooth =>
        exact (plain_genuine_critical_iff hf.1).mpr
          (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr (Or.inr hf.2))))))

/-! ### Per-family type lemmas -/

theorem plain_zeroMeasure_global {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .zeroMeasure β s c θ) : IsGlobalMin .plainRelu β s c θ :=
  isGlobalMin_of_L_eq_zero (L_eq_zero_of_zeroMeasure (P_eq_zero_of_plainZeroField h.1)
    (selfEnergy_eq_zero_of_plainZeroField h.1) h.2)

theorem plain_oneRayExact_global {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .oneRayExact β s c θ) : IsGlobalMin .plainRelu β s c θ := by
  sorry

/-- `stp:census-plain-reduction`: the antipodal cancellation `(a, −a)` at the
opposite ray is not a local minimum. -/
theorem plain_oneRayAntipodalCancel_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .oneRayAntipodalCancel β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

/-- `stp:census-plain-antipodal-collision`: interior collisions are saddles. -/
theorem plain_antipodalInteriorCollision_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .antipodalInteriorCollision β s c θ) :
    IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

/-- `stp:census-plain-antipodal-collision`: the endpoint cusp selector `N > 0`
decides minimum versus saddle. -/
theorem plain_antipodalEndpointCollision_type {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .antipodalEndpointCollision β s c θ) :
    plainFamilyType .antipodalEndpointCollision β s c θ := by
  sorry

theorem plain_oppositePerp_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .oppositePerp β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

/-- `stp:census-plain-antipodal-separated`: the symmetric midpoint pair is a saddle. -/
theorem plain_separatedMidpointPair_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .separatedMidpointPair β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

theorem plain_exactFit_global {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .exactFit β s c θ) : IsGlobalMin .plainRelu β s c θ :=
  isGlobalMin_of_L_eq_zero (L_eq_zero_of_exactFit .plainRelu β s c θ h.2)

/-- `stp:census-plain-gap-collision`: the collision selector decides. -/
theorem plain_collisionRay_type {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .collisionRay β s c θ) : plainFamilyType .collisionRay β s c θ := by
  sorry

/-- `stp:census-plain-gap-collision`: the opposite bisector is a saddle. -/
theorem plain_oppositeBisector_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .oppositeBisector β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

theorem plain_oneDead_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .oneDead β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  have hcrit : IsCritical .plainRelu β s c θ :=
    (plain_two_student_classification β s c θ).mpr ⟨.oneDead, h⟩
  rw [isTopologicalSaddle_iff_not_isLocalMin hcrit]
  obtain ⟨i, j, _, hi, _⟩ := h.2
  exact dead_not_isLocalMin .plainRelu β s c θ (Or.inr ⟨rfl, h.1⟩) ⟨i, hi⟩

theorem plain_bothDead_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .bothDead β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  have hcrit : IsCritical .plainRelu β s c θ :=
    (plain_two_student_classification β s c θ).mpr ⟨.bothDead, h⟩
  rw [isTopologicalSaddle_iff_not_isLocalMin hcrit]
  exact dead_not_isLocalMin .plainRelu β s c θ (Or.inr ⟨rfl, h.1⟩) ⟨0, h.2.1⟩

/-- `clm:census-plain-ghost-diagonal`: at a ghost corner the diagonal entry of
the raw angular block is negative; the configuration is a saddle. -/
theorem plain_separatedGhost_saddle {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .separatedGhost β s c θ) : IsTopologicalSaddle .plainRelu β s c θ := by
  sorry

/-- `stp:census-plain-gap-separated` (Schur clause), `thm-separated-index-law`:
the Schur complement types the smooth separated pair; positive definite is a
strict (spurious) minimum, a negative direction is a saddle, and the rank-one
boundary is decided one order up. -/
theorem plain_separatedSmooth_type {β s c θ : Fin 2 → ℝ}
    (h : PlainMember .separatedSmooth β s c θ) : SchurType β s c θ := by
  sorry

/-- **Each plain-ReLU row carries its displayed label**
(`thm:plain-two-teacher-census`, type clause). -/
theorem plain_two_student_type (f : PlainFamily) (β s c θ : Fin 2 → ℝ)
    (h : PlainMember f β s c θ) : plainFamilyType f β s c θ := by
  cases f with
  | zeroMeasure => exact plain_zeroMeasure_global h
  | oneRayExact => exact plain_oneRayExact_global h
  | oneRayAntipodalCancel => exact plain_oneRayAntipodalCancel_saddle h
  | antipodalInteriorCollision => exact plain_antipodalInteriorCollision_saddle h
  | antipodalEndpointCollision => exact plain_antipodalEndpointCollision_type h
  | oppositePerp => exact plain_oppositePerp_saddle h
  | separatedMidpointPair => exact plain_separatedMidpointPair_saddle h
  | exactFit => exact plain_exactFit_global h
  | collisionRay => exact plain_collisionRay_type h
  | oppositeBisector => exact plain_oppositeBisector_saddle h
  | oneDead => exact plain_oneDead_saddle h
  | bothDead => exact plain_bothDead_saddle h
  | separatedGhost => exact plain_separatedGhost_saddle h
  | separatedSmooth => exact plain_separatedSmooth_type h

/-- Every plain-ReLU critical point is a global minimum, a spurious local
minimum, or a topological saddle. -/
theorem plain_variationalType {β s c θ : Fin 2 → ℝ} (h : IsCriticalL_R β s c θ) :
    VariationalType .plainRelu β s c θ :=
  variationalType_of_isCritical h

/-- `clm:census-plain-load-sign`: the balance scalar `𝖦(D,·)` is negative on
`(0, D)`, positive on `(D, 2π)`, and vanishes in `[0, 2π)` exactly at `0, D`
(with `𝖦(D,x) = (π Φ_R(x−D) − Φ_R(D) Φ_R(x)) H_R(D) − Δ_R(D) H_R(x)`,
`eq-separated-scalar`). -/
theorem plain_load_sign {D : ℝ} (hD0 : 0 < D) (hDπ : D < 2 * π) (hDne : D ≠ π) :
    (∀ x, 0 < x → x < D →
      (π * phiCosJ (x - D) - phiCosJ D * phiCosJ x) * couplingHJ D -
        Δ .plainRelu D * couplingHJ x < 0) ∧
    (∀ x, D < x → x < 2 * π →
      0 < (π * phiCosJ (x - D) - phiCosJ D * phiCosJ x) * couplingHJ D -
        Δ .plainRelu D * couplingHJ x) := by
  sorry

end PlainFamilies

end TwoTeacherCensus
