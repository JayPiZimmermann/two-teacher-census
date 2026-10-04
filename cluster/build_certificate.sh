#!/usr/bin/env bash
# Build the whole development, audit every theorem's axioms, and write an
# immutable record of what was checked into BUILD_CERTIFICATE.txt.
# Modeled on the paper repository's scripts/build_certificate.sh.
#
# The certificate is a convenience record, not a proof: the trusted evidence is
# the Lean source together with the Lean kernel.  It states which commit, which
# toolchain and which dependency lockfile were used, and whether `lake build`
# and the axiom audit succeeded, and how many `sorry` remain.
#
# Usage:  ./cluster/build_certificate.sh   (from the repository root)

set -uo pipefail
cd "$(dirname "$0")/.."
OUT=BUILD_CERTIFICATE.txt

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1
  else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

if ! git rev-parse --verify -q HEAD^{commit} >/dev/null 2>&1; then
  echo "error: this repository has no commit yet; commit first." >&2
  exit 2
fi

commit="$(git rev-parse HEAD)"
tree_state="clean"
if [ -n "$(git status --porcelain 2>/dev/null | grep -v " ${OUT}$")" ]; then
  tree_state="dirty (uncommitted changes present)"
fi

lean_version="$(lean --version 2>/dev/null | head -1)"
lake_version="$(lake --version 2>/dev/null | head -1)"
toolchain="$(cat lean-toolchain)"
manifest_hash="$(sha256 lake-manifest.json)"
lakefile_hash="$(sha256 lakefile.lean)"

echo "== lake build"
build_log="$(mktemp)"
lake build 2>&1 | tee "$build_log"
build_status=${PIPESTATUS[0]}
# `lake build` succeeds on a file that still contains `sorry`; Lean only warns.
sorry_count="$(grep -c "declaration uses 'sorry'" "$build_log" || true)"
modules="$(ls TwoTeacherCensus/*.lean | wc -l)"

echo "== axioms of every theorem"
sweep=cluster/AxiomSweep.lean
cluster/gen_axiom_sweep.sh > "$sweep"
statements="$(grep -c '^#print axioms' "$sweep")"
sweep_out="$(lake env lean "./$sweep" 2>&1)"
rm -f "$sweep"
reported="$(printf '%s\n' "$sweep_out" | grep -c 'depends on axioms\|does not depend on any axioms')"
used="$(printf '%s ' "$sweep_out" | grep -oE '\[[^]]*\]' | tr -d ' []' | tr ',' '\n' | sort -u | paste -sd' ')"
unexpected="$(printf '%s\n' "$sweep_out" | tr '\n' ' ' | grep -oE '\[[^]]*\]' | tr -d ' []' | tr ',' '\n' \
  | grep -vE '^(propext|Classical\.choice|Quot\.sound)$' | sort -u | paste -sd' ')"
sorry_dependent="$(printf '%s\n' "$sweep_out" | tr '\n' ' ' | grep -oE "'[^']+' depends on axioms: \[[^]]*sorryAx[^]]*\]" | wc -l)"

{
  echo "BUILD CERTIFICATE -- two-teacher-census"
  echo
  echo "This file records one complete check of the development.  It is not a"
  echo "cryptographic proof of correctness: the trusted evidence is the Lean"
  echo "source in this repository together with the Lean kernel.  What it fixes"
  echo "is exactly which source, toolchain and dependencies were checked."
  echo
  echo "date (UTC)             : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "commit                 : ${commit}"
  echo "working tree           : ${tree_state}"
  echo "toolchain (lean-toolchain): ${toolchain}"
  echo "lean --version         : ${lean_version}"
  echo "lake --version         : ${lake_version}"
  echo "lake-manifest.json     : sha256 ${manifest_hash}"
  echo "lakefile.lean          : sha256 ${lakefile_hash}"
  echo
  echo "lake build             : exit ${build_status}$( [ "$build_status" -eq 0 ] && echo '  (success)' )"
  echo "modules compiled       : ${modules} (TwoTeacherCensus/*.lean; TwoTeacherCensus.lean imports every one)"
  echo "'sorry' warnings       : ${sorry_count}"
  echo "theorems audited       : ${reported} of ${statements}"
  echo "theorems using sorryAx : ${sorry_dependent}"
  echo "axioms used            : ${used}"
  echo "unexpected axioms      : ${unexpected:-none}"
  echo
  echo "The axiom sweep prints the transitive axiom dependencies of every theorem"
  echo "and lemma; anything beyond propext, Classical.choice and Quot.sound -- in"
  echo "particular sorryAx and Lean.ofReduceBool -- is listed as unexpected."
  echo
  echo "Reproduce with:"
  echo "    ./cluster/build_certificate.sh"
} > "$OUT"

echo
cat "$OUT"
rm -f "$build_log"
[ "$build_status" -eq 0 ] && [ "$sorry_count" -eq 0 ] && [ -z "$unexpected" ] \
  && [ "$reported" -eq "$statements" ]
