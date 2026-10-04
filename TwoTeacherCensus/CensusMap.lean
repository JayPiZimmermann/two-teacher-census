import TwoTeacherCensus.Classification

/-!
# The census map

The census of a two-atom teacher as a function of the teacher: the map
coordinates `(β̂, y)` of `eq-map-mass-coordinate`, the three mass-free
determinants `𝒲_τ, 𝒲_P, 𝒲_W` of `eq-map-wronskians`, the sign bridge
`eq-map-sign-bridge`, and the centered census-map theorem
`thm-map-centered-complete` (nine boundary curves, ten faces, four census
values, five positive-area regions).  The plain-ReLU map is stated as a
conjecture with its certified parts named (`certificates/COMPLETENESS.md`).

The census is read with multiplicity on the colliding rows: a `CensusValue`
records the set of entries and the number of distinct collision directions
carrying a spurious minimum of each split sign, so that `fit + coin⁻` and
`fit + 2 coin⁻` are different values.
-/

noncomputable section

open Real Set Filter Topology Classical
open scoped BigOperators

namespace TwoTeacherCensus

/-! ## Map coordinates (`eq-map-mass-coordinate`) -/

/-- The teacher directions `(0, β̂)` of a map point. -/
def mapTeacher (b : ℝ) : Fin 2 → ℝ := ![0, b]

/-- The mass pair `(sin ψ, cos ψ)`, `ψ = (y+1)π/2`, of a map point. -/
def mapMass (y : ℝ) : Fin 2 → ℝ := ![sin ((y + 1) * π / 2), cos ((y + 1) * π / 2)]

/-- `y ∈ (−1, 0)`: both masses positive (after the global sign flip, every
same-sign teacher). -/
theorem mapMass_pos_of_neg {y : ℝ} (h1 : -1 < y) (h0 : y < 0) :
    0 < mapMass y 0 ∧ 0 < mapMass y 1 := by
  unfold mapMass
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
  have hψ0 : 0 < (y + 1) * π / 2 := by
    have := Real.pi_pos; nlinarith
  have hψ1 : (y + 1) * π / 2 < π / 2 := by
    have := Real.pi_pos; nlinarith
  exact ⟨sin_pos_of_pos_of_lt_pi hψ0 (by linarith),
    cos_pos_of_mem_Ioo ⟨by linarith, hψ1⟩⟩

/-- `y ∈ (0, 1)`: masses of opposite signs. -/
theorem mapMass_mixed_of_pos {y : ℝ} (h0 : 0 < y) (h1 : y < 1) :
    0 < mapMass y 0 ∧ mapMass y 1 < 0 := by
  unfold mapMass
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons]
  have hψ0 : π / 2 < (y + 1) * π / 2 := by
    have := Real.pi_pos; nlinarith
  have hψ1 : (y + 1) * π / 2 < π := by
    have := Real.pi_pos; nlinarith
  exact ⟨sin_pos_of_pos_of_lt_pi (by linarith [Real.pi_pos]) hψ1,
    cos_neg_of_pi_div_two_lt_of_lt hψ0 (by linarith [Real.pi_pos])⟩

/-! ## The three mass-free determinants (`eq-map-wronskians`) -/

/-- `𝒲_τ(β,t) = H_q(t) H_q'(t−β) − H_q(t−β) H_q'(t)`. -/
def wTau (q : Model) (b t : ℝ) : ℝ :=
  Hq q t * Hslope q (t - b) - Hq q (t - b) * Hslope q t

/-- `𝒲_P(β,t) = Φ_q(t) H_q(t−β) − Φ_q(t−β) H_q(t)`. -/
def wPot (q : Model) (b t : ℝ) : ℝ :=
  Φ q t * Hq q (t - b) - Φ q (t - b) * Hq q t

/-- `𝒲_W(β,t) = |sin t| H_q(t−β) − |sin(t−β)| H_q(t)`. -/
def wWeight (q : Model) (b t : ℝ) : ℝ :=
  |sin t| * Hq q (t - b) - |sin (t - b)| * Hq q t

/-- `eq-map-wronskian-relation`: `𝒲_τ + 𝒲_P = 2 𝒲_W`, from the kernel equation
`H_q' = Φ_q − 2|sin|`. -/
theorem wTau_add_wPot (q : Model) (b t : ℝ) :
    wTau q b t + wPot q b t = 2 * wWeight q b t := by
  unfold wTau wPot wWeight Hslope
  ring

/-- `eq-map-mass-line`: a torque root at which the two coupling atoms do not
both vanish determines the mass direction up to scale. -/
theorem massLine_of_torqueRoot (q : Model) {b s0 s1 t : ℝ}
    (hA : s0 * Hq q t + s1 * Hq q (t - b) = 0)
    (hne : Hq q t ≠ 0 ∨ Hq q (t - b) ≠ 0) :
    ∃ l : ℝ, s0 = l * (-Hq q (t - b)) ∧ s1 = l * Hq q t := by
  rcases hne with h | h
  · refine ⟨s1 / Hq q t, ?_, ?_⟩
    · field_simp; linarith
    · field_simp
  · refine ⟨-s0 / Hq q (t - b), ?_, ?_⟩
    · field_simp
    · field_simp
      have : s1 * Hq q (t - b) = -(s0 * Hq q t) := by linarith
      linear_combination this

/-- `eq-map-sign-bridge`: on the mass line `(s_0,s_1) = λ(−H_q(t−β), H_q(t))`,
`P_q(t) = −λ 𝒲_P`, `W(t) = −λ 𝒲_W`, `τ_q(t) = λ 𝒲_τ`. -/
theorem signBridge (q : Model) (b l t : ℝ) :
    P q (mapTeacher b) ![l * (-Hq q (t - b)), l * Hq q t] t = -l * wPot q b t ∧
    W (mapTeacher b) ![l * (-Hq q (t - b)), l * Hq q t] t = -l * wWeight q b t ∧
    τ q (mapTeacher b) ![l * (-Hq q (t - b)), l * Hq q t] t = l * wTau q b t := by
  have hrel := wTau_add_wPot q b t
  unfold τ P W wPot wWeight mapTeacher
  simp only [Fin.sum_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    sub_zero]
  refine ⟨by ring, by ring, ?_⟩
  unfold wPot wWeight at hrel
  linear_combination (-l) * hrel

/-! ## Census values with multiplicity -/

/-- Collision directions in one period carrying a spurious minimum of the given
split type. -/
def collisionTrapDirections (q : Model) (β s : Fin 2 → ℝ) (m : MinType) : Set ℝ :=
  {t | 0 ≤ t ∧ t < period q ∧
    ∃ c : Fin 2 → ℝ, IsLocalMinL q β s c ![t, t] ∧ censusEntry q β s c ![t, t] = (.colliding, m)}

/-- A census value: the set of entries and the colliding multiplicities. -/
structure CensusValue where
  entries : Set (FamilyKind × MinType)
  coinPositive : ℕ
  coinMixed : ℕ

/-- The census value of a teacher. -/
def censusValue (q : Model) (β s : Fin 2 → ℝ) : CensusValue :=
  ⟨Census q β s,
    (collisionTrapDirections q β s .trapPositive).ncard,
    (collisionTrapDirections q β s .trapMixed).ncard⟩

/-- `fit + coin⁻`. -/
def fitCoinMinus : CensusValue := ⟨{(.fit, .global), (.colliding, .trapMixed)}, 0, 1⟩
/-- `fit + 2 coin⁻`. -/
def fitTwoCoinMinus : CensusValue := ⟨{(.fit, .global), (.colliding, .trapMixed)}, 0, 2⟩
/-- `fit + coin⁺`. -/
def fitCoinPlus : CensusValue := ⟨{(.fit, .global), (.colliding, .trapPositive)}, 1, 0⟩
/-- `fit + 2 coin⁺`. -/
def fitTwoCoinPlus : CensusValue := ⟨{(.fit, .global), (.colliding, .trapPositive)}, 2, 0⟩

/-! ## Boundary curves of the centered map -/

/-- Torque roots of a map teacher in one period. -/
def torqueRoots (q : Model) (b y : ℝ) : Set ℝ :=
  {t | 0 ≤ t ∧ t < period q ∧ A q (mapTeacher b) (mapMass y) t = 0}

/-- The lens: map points with a degenerate torque root (`𝒲_τ = 0` at a root). -/
def LensCurve : Set (ℝ × ℝ) :=
  {p | ∃ t ∈ torqueRoots .centered p.1 p.2, wTau .centered p.1 t = 0}

/-- The potential curve: map points with a torque root at which `P_C` vanishes. -/
def PotentialCurve : Set (ℝ × ℝ) :=
  {p | ∃ t ∈ torqueRoots .centered p.1 p.2, wPot .centered p.1 t = 0}

/-- The nine boundary curves of the centered map: lens, potential curve, the
lines `y = 0`, `y = ±1` and the columns `β ∈ {0, π/2, π}`. -/
def CenteredBoundary : Set (ℝ × ℝ) :=
  LensCurve ∪ PotentialCurve ∪ {p | p.2 = 0 ∨ p.2 = 1 ∨ p.2 = -1} ∪
    {p | p.1 = 0 ∨ p.1 = π / 2 ∨ p.1 = π}

/-- `prop-map-boundary-elimination`: at a torque root of a nonzero teacher, a
degenerate root forces `𝒲_τ = 0`, a vanishing potential forces `𝒲_P = 0`, and a
vanishing weight forces `𝒲_W = 0`. -/
theorem boundary_elimination (q : Model) {b s0 s1 t : ℝ} (hs : s0 ≠ 0 ∨ s1 ≠ 0)
    (hA : s0 * Hq q t + s1 * Hq q (t - b) = 0) :
    (s0 * Hslope q t + s1 * Hslope q (t - b) = 0 → wTau q b t = 0) ∧
    (s0 * Φ q t + s1 * Φ q (t - b) = 0 → wPot q b t = 0) ∧
    (s0 * |sin t| + s1 * |sin (t - b)| = 0 → wWeight q b t = 0) := by
  have key : ∀ u v : ℝ, s0 * u + s1 * v = 0 →
      u * Hq q (t - b) - v * Hq q t = 0 := by
    intro u v huv
    rcases hs with h | h
    · have := congrArg (fun z => z * s0) (show (u * Hq q (t - b) - v * Hq q t) * s0 =
        (s0 * u + s1 * v) * Hq q (t - b) - (s0 * Hq q t + s1 * Hq q (t - b)) * v by ring)
      have h2 : (u * Hq q (t - b) - v * Hq q t) * s0 = 0 := by
        rw [show (u * Hq q (t - b) - v * Hq q t) * s0 =
          (s0 * u + s1 * v) * Hq q (t - b) - (s0 * Hq q t + s1 * Hq q (t - b)) * v by ring,
          huv, hA]; ring
      exact (mul_eq_zero.mp h2).resolve_right h
    · have h2 : (u * Hq q (t - b) - v * Hq q t) * s1 = 0 := by
        rw [show (u * Hq q (t - b) - v * Hq q t) * s1 =
          (s0 * Hq q t + s1 * Hq q (t - b)) * u - (s0 * u + s1 * v) * Hq q t by ring,
          huv, hA]; ring
      exact (mul_eq_zero.mp h2).resolve_right h
  refine ⟨fun h => ?_, fun h => ?_, fun h => ?_⟩
  · have := key _ _ h; unfold wTau; linarith
  · have := key _ _ h; unfold wPot; linarith
  · have := key _ _ h; unfold wWeight; linarith

/-! ## Root counts (`thm-map-centered-counts`) -/

/-- The constant `β*_{C,2} ∈ (π/2, π)`: the root of `(x/2) tan(x/2) = 1`. -/
def IsBetaStarC2 (x : ℝ) : Prop := π / 2 < x ∧ x < π ∧ (x / 2) * tan (x / 2) = 1

/-- `eq-map-centered-betastar`: the constant exists uniquely. -/
theorem existsUnique_betaStarC2 : ∃! x, IsBetaStarC2 x := by
  sorry

/-- `β*_{C,2}`. -/
def βstarC2 : ℝ := Classical.choose existsUnique_betaStarC2.exists

/-- `β*_{C,1} = π − β*_{C,2}`, equivalently `2u*` with `tan u* = π/2 − u*`. -/
def βstarC1 : ℝ := π - βstarC2

/-- `thm-map-centered-counts`: per period, `𝒲_W` vanishes exactly at `t = 0, β`;
`𝒲_P` has exactly two roots; `𝒲_τ` has no root for
`β ∉ [β*_{C,1}, β*_{C,2}]`, one double root at the endpoints, two roots inside. -/
theorem centered_root_counts {b : ℝ} (hb0 : 0 < b) (hbπ : b < π) :
    ({t | 0 ≤ t ∧ t < π ∧ wWeight .centered b t = 0} = {0, b}) ∧
    ({t | 0 ≤ t ∧ t < π ∧ wPot .centered b t = 0}.ncard = 2) ∧
    ((b < βstarC1 ∨ βstarC2 < b) → {t | 0 ≤ t ∧ t < π ∧ wTau .centered b t = 0} = ∅) ∧
    ((b = βstarC1 ∨ b = βstarC2) → {t | 0 ≤ t ∧ t < π ∧ wTau .centered b t = 0}.ncard = 1) ∧
    ((βstarC1 < b ∧ b < βstarC2) → {t | 0 ≤ t ∧ t < π ∧ wTau .centered b t = 0}.ncard = 2) := by
  sorry

/-! ## The faces -/

/-- Lower and upper branch of the lens over the window `[β*_{C,1}, β*_{C,2}]`
(the lens lies in the sector `y ∈ (−1, 0)` and is symmetric about `y = −½`). -/
def lensLower (b : ℝ) : ℝ := sInf {y | -1 < y ∧ y < 0 ∧ (b, y) ∈ LensCurve}
def lensUpper (b : ℝ) : ℝ := sSup {y | -1 < y ∧ y < 0 ∧ (b, y) ∈ LensCurve}

/-- Lower and upper branch of the potential curve over `(0, π)` (in the sector
`y ∈ (0, 1)`, symmetric about `y = ½`). -/
def potLower (b : ℝ) : ℝ := sInf {y | 0 < y ∧ y < 1 ∧ (b, y) ∈ PotentialCurve}
def potUpper (b : ℝ) : ℝ := sSup {y | 0 < y ∧ y < 1 ∧ (b, y) ∈ PotentialCurve}

/-- Strictly inside the lens. -/
def InsideLens (b y : ℝ) : Prop :=
  βstarC1 < b ∧ b < βstarC2 ∧ lensLower b < y ∧ y < lensUpper b

/-- Strictly outside the closed lens. -/
def OutsideLens (b y : ℝ) : Prop :=
  ¬ (βstarC1 ≤ b ∧ b ≤ βstarC2 ∧ lensLower b ≤ y ∧ y ≤ lensUpper b)

/-- The ten faces of the centered map (`thm-map-centered-complete`). -/
inductive CenteredFace where
  | F1 | F2 | F3 | F4 | F5 | F6 | F7 | F8 | F9 | F10
  deriving DecidableEq

/-- Face membership, as the explicit open regions of the table. -/
def CenteredFaceMem : CenteredFace → ℝ → ℝ → Prop
  | .F1, b, y => 0 < b ∧ b < π / 2 ∧ -1 < y ∧ y < 0 ∧ OutsideLens b y
  | .F2, b, y => 0 < b ∧ b < π / 2 ∧ InsideLens b y
  | .F3, b, y => π / 2 < b ∧ b < π ∧ InsideLens b y
  | .F4, b, y => π / 2 < b ∧ b < π ∧ -1 < y ∧ y < 0 ∧ OutsideLens b y
  | .F5, b, y => 0 < b ∧ b < π / 2 ∧ 0 < y ∧ y < potLower b
  | .F6, b, y => 0 < b ∧ b < π / 2 ∧ potLower b < y ∧ y < potUpper b
  | .F7, b, y => 0 < b ∧ b < π / 2 ∧ potUpper b < y ∧ y < 1
  | .F8, b, y => π / 2 < b ∧ b < π ∧ 0 < y ∧ y < potLower b
  | .F9, b, y => π / 2 < b ∧ b < π ∧ potLower b < y ∧ y < potUpper b
  | .F10, b, y => π / 2 < b ∧ b < π ∧ potUpper b < y ∧ y < 1

/-- The census of each face (table of `thm-map-centered-complete`). -/
def centeredFaceCensus : CenteredFace → CensusValue
  | .F1 => fitCoinMinus
  | .F2 => fitTwoCoinMinus
  | .F3 => fitTwoCoinMinus
  | .F4 => fitCoinMinus
  | .F5 => fitCoinPlus
  | .F6 => fitTwoCoinPlus
  | .F7 => fitCoinPlus
  | .F8 => fitCoinPlus
  | .F9 => fitTwoCoinPlus
  | .F10 => fitCoinPlus

/-- The number of torque roots of each face (`n` in the table). -/
def centeredFaceRootCount : CenteredFace → ℕ
  | .F2 => 4
  | .F3 => 4
  | _ => 2

/-- Structure of the two curves: the lens is a closed curve over
`[β*_{C,1}, β*_{C,2}]` inside the sector `y ∈ (−1, 0)` with two branches meeting
only at the endpoints on `y = −½`; the potential curve has two disjoint branches
over all of `(0, π)` inside the sector `y ∈ (0, 1)`.  (`appendix-map-centered`,
`eq-map-potential-small-gap`.) -/
theorem centered_curves_structure :
    (∀ b, βstarC1 < b → b < βstarC2 → -1 < lensLower b ∧ lensLower b < -1 / 2 ∧
      -1 / 2 < lensUpper b ∧ lensUpper b < 0 ∧
      (∀ y, -1 < y → y < 0 → ((b, y) ∈ LensCurve ↔ y = lensLower b ∨ y = lensUpper b))) ∧
    (∀ b, 0 < b → b < π → 0 < potLower b ∧ potLower b < 1 / 2 ∧ 1 / 2 < potUpper b ∧
      potUpper b < 1 ∧
      (∀ y, 0 < y → y < 1 → ((b, y) ∈ PotentialCurve ↔ y = potLower b ∨ y = potUpper b))) ∧
    (∀ b y, -1 < y → y < 0 → (b, y) ∈ LensCurve → βstarC1 ≤ b ∧ b ≤ βstarC2) := by
  sorry

/-- **The centered census map, exactly** (`thm-map-centered-complete`): the
census is constant on each of the ten faces with the tabulated value, and the
number of torque roots is as tabulated. -/
theorem centered_census_map (f : CenteredFace) {b y : ℝ} (h : CenteredFaceMem f b y) :
    censusValue .centered (mapTeacher b) (mapMass y) = centeredFaceCensus f ∧
    (torqueRoots .centered b y).ncard = centeredFaceRootCount f := by
  sorry

/-- The census changes only across the nine curves: off `CenteredBoundary` it is
locally constant (`thm-param-zero-count`, `appendix-map-boundaries`). -/
theorem centered_census_locally_constant {b y : ℝ} (h : (b, y) ∉ CenteredBoundary)
    (hb : 0 < b ∧ b < π) (hy : -1 < y ∧ y < 1) :
    ∀ᶠ p : ℝ × ℝ in 𝓝 (b, y),
      censusValue .centered (mapTeacher p.1) (mapMass p.2) =
        censusValue .centered (mapTeacher b) (mapMass y) := by
  sorry

/-- Four distinct census values occur, and merging across the measure-zero
column `β = π/2` leaves five regions: `F1F4, F2F3, F5F8, F6F9, F7F10`. -/
theorem centered_census_four_values :
    (∀ f, centeredFaceCensus f = fitCoinMinus ∨ centeredFaceCensus f = fitTwoCoinMinus ∨
      centeredFaceCensus f = fitCoinPlus ∨ centeredFaceCensus f = fitTwoCoinPlus) ∧
    centeredFaceCensus .F1 = centeredFaceCensus .F4 ∧
    centeredFaceCensus .F2 = centeredFaceCensus .F3 ∧
    centeredFaceCensus .F5 = centeredFaceCensus .F8 ∧
    centeredFaceCensus .F6 = centeredFaceCensus .F9 ∧
    centeredFaceCensus .F7 = centeredFaceCensus .F10 := by
  refine ⟨fun f => ?_, rfl, rfl, rfl, rfl, rfl⟩
  cases f <;> simp [centeredFaceCensus]

/-! ## The plain-ReLU map (partially certified) -/

/-- Plain-ReLU boundary elimination (`prop-map-boundary-elimination-plain`) is
the `q = .plainRelu` instance of `boundary_elimination`. -/
theorem plain_boundary_elimination {b s0 s1 t : ℝ} (hs : s0 ≠ 0 ∨ s1 ≠ 0)
    (hA : s0 * Hq .plainRelu t + s1 * Hq .plainRelu (t - b) = 0) :
    (s0 * Hslope .plainRelu t + s1 * Hslope .plainRelu (t - b) = 0 → wTau .plainRelu b t = 0) ∧
    (s0 * Φ .plainRelu t + s1 * Φ .plainRelu (t - b) = 0 → wPot .plainRelu b t = 0) ∧
    (s0 * |sin t| + s1 * |sin (t - b)| = 0 → wWeight .plainRelu b t = 0) :=
  boundary_elimination .plainRelu hs hA

/-- The constant `β*_{R,1} = 2u*_R`, `tan u*_R = π − u*_R`, `u*_R ∈ (0, π/2)`. -/
def IsUStarR (u : ℝ) : Prop := 0 < u ∧ u < π / 2 ∧ tan u = π - u

theorem existsUnique_uStarR : ∃! u, IsUStarR u := by
  sorry

def βstarR1 : ℝ := 2 * Classical.choose existsUnique_uStarR.exists

/-- `thm-map-plain-relu-counts` (certified in `certificates/census_map_noncentered/
analytic_counts.json` for every `β ∈ (0, 2π) \ {π}`; proved here only as a
statement): `𝒲_P` has no interior roots, its only zeros per period being
`t = π, β + π`; `𝒲_τ` has interior roots only on `(0, β)`, none for
`β < β*_{R,1}`, a double root at `β/2` for `β = β*_{R,1}`, and at least one on
each side of the midpoint for `β > β*_{R,1}`. -/
theorem plain_root_counts {b : ℝ} (hb0 : 0 < b) (hbπ : b < π) :
    ({t | 0 ≤ t ∧ t < 2 * π ∧ wPot .plainRelu b t = 0} = {π, b + π}) ∧
    (∀ t, 0 ≤ t → t < 2 * π → wTau .plainRelu b t = 0 → t = π ∨ t = b + π ∨ (0 < t ∧ t < b)) ∧
    (b < βstarR1 → ∀ t, 0 < t → t < b → wTau .plainRelu b t ≠ 0) ∧
    (b = βstarR1 → ∀ t, 0 < t → t < b → (wTau .plainRelu b t = 0 ↔ t = b / 2)) ∧
    (βstarR1 < b → (∃ t, 0 < t ∧ t < b / 2 ∧ wTau .plainRelu b t = 0) ∧
      (∃ t, b / 2 < t ∧ t < b ∧ wTau .plainRelu b t = 0)) := by
  sorry

/-- The colliding half of the plain-ReLU census is a function of the sign chart
(`eq-map-selector-signs`): off the degenerate locus `𝒲_τ = 0`, a torque root `t`
carries a spurious collision minimum iff `𝒲_τ 𝒲_P < 0` there, with split sign
`sign(c_0 c_1) = sign(𝒲_τ 𝒲_W)`.  Follows from `plain_collisionRay_type` and the
sign bridge. -/
theorem plain_colliding_census_of_signChart {b y t : ℝ} (hb : 0 < b ∧ b < π)
    (hy : -1 < y ∧ y < 1 ∧ y ≠ 0) (ht : t ∈ torqueRoots .plainRelu b y)
    (hτ : wTau .plainRelu b t ≠ 0) :
    ((∃ c : Fin 2 → ℝ, IsSpuriousMin .plainRelu (mapTeacher b) (mapMass y) c ![t, t]) ↔
      wTau .plainRelu b t * wPot .plainRelu b t < 0) := by
  sorry

/-- **Conjecture (plain-ReLU census map).**  The plain-ReLU census is locally
constant off the union of: the torque lens (`𝒲_τ = 0` at a torque root), the
lines `y ∈ {0, ±1}`, the columns `β ∈ {0, π, 2π}`, and the separated type
boundary (`eq-map-separated-type-boundary`, the index-degeneracy locus of the
balance map), and the drawn seven-label arrangement of
`certificates/census_map_noncentered/faces_fixed.json` is its face structure.
Status (`certificates/COMPLETENESS.md`): coincident root counts and the colliding
half are certified for every gap; the separated half is certified only at the
22 exhaustively enumerated teachers, on the F4 rectangle `[0.7,0.8]×[−0.45,−0.35]`
and on the 68-tile atlas `[0.5,2.2]×[−0.45,−0.1]`; the F5 label is known NOT to
carry one constant family count.  This statement is a conjecture and is not
expected to be closed by the cluster run. -/
theorem plain_census_map_conjecture :
    ∃ B : Set (ℝ × ℝ), (LensCurve ∪ {p | p.2 = 0 ∨ p.2 = 1 ∨ p.2 = -1} ∪
        {p | p.1 = 0 ∨ p.1 = π ∨ p.1 = 2 * π}) ⊆ B ∧
      ∀ b y, 0 < b → b < 2 * π → -1 < y → y < 1 → (b, y) ∉ B →
        ∀ᶠ p : ℝ × ℝ in 𝓝 (b, y),
          censusValue .plainRelu (mapTeacher p.1) (mapMass p.2) =
            censusValue .plainRelu (mapTeacher b) (mapMass y) := by
  sorry

end TwoTeacherCensus
