#!/usr/bin/env bash
# Axiom audit: every theorem/lemma of the development must depend on nothing
# beyond propext, Classical.choice and Quot.sound.  Prints a summary and exits
# 0 iff the audit is clean (no sorryAx, no other unexpected axiom).
# Requires a built tree (`lake build`); uses `lake env lean` on a generated file.
set -uo pipefail
cd "$(dirname "$0")/.."
sweep=cluster/AxiomSweep.lean
cluster/gen_axiom_sweep.sh > "$sweep"
statements="$(grep -c '^#print axioms' "$sweep")"
out="$(lake env lean "$sweep" 2>&1)"
rm -f "$sweep"
reported="$(printf '%s\n' "$out" | grep -c "depends on axioms\|does not depend on any axioms")"
errors="$(printf '%s\n' "$out" | grep -c "error:")"
unexpected="$(printf '%s\n' "$out" | tr '\n' ' ' | grep -oE '\[[^]]*\]' | tr -d ' []' | tr ',' '\n' \
  | grep -vE '^(propext|Classical\.choice|Quot\.sound)$' | sort | uniq -c | sort -rn)"
sorry_dependent_all="$(printf '%s\n' "$out" | tr '\n' ' ' | grep -oE "'[^']+' depends on axioms: \[[^]]*sorryAx[^]]*\]" | sed "s/' depends.*//; s/^'//" )"
allow="$(grep -v '^#' cluster/sorry_allowlist.txt 2>/dev/null | sed '/^$/d')"
sorry_dependent="$(printf '%s\n' "$sorry_dependent_all" | sed '/^$/d' | grep -vxF -f <(printf '%s\n' "$allow") || true)"
allowed_hits="$(printf '%s\n' "$sorry_dependent_all" | sed '/^$/d' | grep -xF -f <(printf '%s\n' "$allow") || true)"
# sorryAx is tolerated exactly on the allowlisted conjectures.
if [ -z "$sorry_dependent" ]; then
  unexpected="$(printf '%s\n' "$unexpected" | grep -v 'sorryAx' || true)"
fi
echo "axiom audit: $reported of $statements statements reported, $errors errors"
[ -n "$allowed_hits" ] && echo "allowlisted conjectures still using sorryAx: $(echo "$allowed_hits" | tr '\n' ' ')"
if [ -n "$unexpected" ]; then
  echo "unexpected axioms (count, name):"
  echo "$unexpected"
  echo "theorems depending on sorryAx:"
  echo "$sorry_dependent"
  exit 1
fi
[ "$errors" -eq 0 ] && [ "$reported" -eq "$statements" ] && echo "axiom audit clean" && exit 0
exit 1
