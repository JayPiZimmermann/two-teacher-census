#!/bin/bash
# Pack and REPLAY every emitted certificate that has not been replayed yet.
# Prints one line per witness: leaves / verified / undecided / verdict.
# The batch log is prospecting; only these replay lines are evidence.
cd "$(dirname "$0")"
for f in ivcert_C_*.json; do
  case "$f" in *.bits.json) continue;; esac
  tag="${f#ivcert_}"; tag="${tag%.json}"
  [ -f "replay_${tag}.log" ] && continue
  python3 ivcert_pack.py "$f" > "pack_${tag}.log" 2>&1
  timeout 6000 python3 ivcert_check.py "ivcert_${tag}.bits.json" \
      > "replay_${tag}.log" 2>&1
  echo "EXIT $?" >> "replay_${tag}.log"
done
