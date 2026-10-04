#!/usr/bin/env bash
# Generate cluster/AxiomSweep.lean: one `#print axioms` per theorem/lemma of the
# development.  Usage: cluster/gen_axiom_sweep.sh > cluster/AxiomSweep.lean
set -euo pipefail
cd "$(dirname "$0")/.."
echo "import TwoTeacherCensus"
for f in TwoTeacherCensus/*.lean; do
  grep -oE "^(theorem|lemma) [^ :({]+" "$f" | awk '{print "#print axioms TwoTeacherCensus." $2}'
done
