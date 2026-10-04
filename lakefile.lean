import Lake
open Lake DSL

package «two_teacher_census» where
  -- no package-level options

require mathlib from git
  "https://github.com/leanprover-community/mathlib4.git" @ "v4.4.0"

/-- The planar two-student / two-teacher classification and census map.
`lake build` compiles every module; `TwoTeacherCensus.lean` imports them all. -/
@[default_target]
lean_lib «TwoTeacherCensus» where
  srcDir := "."
  roots := #[`TwoTeacherCensus, `TwoTeacherCensus.Definitions,
    `TwoTeacherCensus.Classification, `TwoTeacherCensus.CensusMap]
  moreLeanArgs := #["-M", "16384"]
