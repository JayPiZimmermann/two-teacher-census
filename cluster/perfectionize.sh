#!/usr/bin/env bash
# Cluster loop: drive a `claude` agent to remove every `sorry` from the Lean
# development, rebuilding and auditing after each pass, committing each
# successful improvement.
#
# Requirements on the cluster: `lake`, `lean` (elan toolchain from
# lean-toolchain), a Mathlib cache (`lake exe cache get`), the `claude` CLI,
# `git`.  Run from the repository root:
#
#     ./cluster/perfectionize.sh
#
# Environment:
#   MAX_ITERS        maximum agent iterations            (default 40)
#   PERMISSION_MODE  claude permission flag               (default "--permission-mode acceptEdits";
#                    set to "--dangerously-skip-permissions" for unattended runs)
#   CLAUDE_BIN       the claude executable                (default "claude")
#   CLAUDE_EXTRA     extra flags for claude -p            (default "")
#   SKIP_CACHE       set to 1 to skip `lake exe cache get`
#   NO_COMMIT        set to 1 to never commit
#   ALLOWED_SORRIES  number of `sorry` warnings tolerated at the end (default 1:
#                    the documented conjecture `plain_census_map_conjecture`,
#                    see cluster/sorry_allowlist.txt)
#
# Stopping condition: `lake build` has 0 errors and no `sorry` warnings beyond
# the allowlisted conjectures (cluster/sorry_allowlist.txt), and the axiom audit
# (cluster/audit_axioms.sh) reports only propext, Classical.choice, Quot.sound
# outside that allowlist.  Each iteration's build log and agent transcript go to cluster/log/.

set -uo pipefail
cd "$(dirname "$0")/.."

MAX_ITERS="${MAX_ITERS:-40}"
PERMISSION_MODE="${PERMISSION_MODE:---permission-mode acceptEdits}"
CLAUDE_BIN="${CLAUDE_BIN:-claude}"
CLAUDE_EXTRA="${CLAUDE_EXTRA:-}"
SKIP_CACHE="${SKIP_CACHE:-0}"
NO_COMMIT="${NO_COMMIT:-0}"
ALLOWED_SORRIES="${ALLOWED_SORRIES:-1}"

mkdir -p cluster/log
stamp() { date -u +%Y%m%dT%H%M%SZ; }
log() { echo "[$(stamp)] $*" | tee -a cluster/log/perfectionize.log; }

count_metrics() {
  # sets ERRORS, SORRIES from cluster/last_build.log
  ERRORS="$(grep -c "error:" cluster/last_build.log || true)"
  SORRIES="$(grep -c "declaration uses 'sorry'" cluster/last_build.log || true)"
  SORRY_LINES="$(grep -c "sorry" TwoTeacherCensus/*.lean | awk -F: '{s+=$2} END {print s+0}')"
}

build_once() {
  lake build 2>&1 | tee cluster/last_build.log >/dev/null
  BUILD_STATUS=${PIPESTATUS[0]}
  count_metrics
}

audit_once() {
  if cluster/audit_axioms.sh > cluster/last_audit.log 2>&1; then AUDIT_OK=1; else AUDIT_OK=0; fi
}

# ---- 1. one-time setup -----------------------------------------------------
if [ "$SKIP_CACHE" != "1" ]; then
  log "lake exe cache get"
  lake exe cache get 2>&1 | tail -3 | tee -a cluster/log/perfectionize.log
fi
log "initial lake build"
build_once
audit_once
log "initial state: build exit $BUILD_STATUS, errors=$ERRORS, sorry warnings=$SORRIES, sorry tokens=$SORRY_LINES, audit_ok=$AUDIT_OK"
cp cluster/last_build.log "cluster/log/build_000.log"

best_sorries="$SORRY_LINES"
best_errors="$ERRORS"

# ---- 2. the loop -----------------------------------------------------------
for ((i = 1; i <= MAX_ITERS; i++)); do
  if [ "$BUILD_STATUS" -eq 0 ] && [ "$ERRORS" -eq 0 ] && [ "$SORRIES" -le "$ALLOWED_SORRIES" ] \
     && [ "$AUDIT_OK" -eq 1 ]; then
    log "DONE: 0 errors, $SORRIES sorry warning(s) (allowed: $ALLOWED_SORRIES, conjectures only), axiom audit clean."
    [ "$NO_COMMIT" = "1" ] || ./cluster/build_certificate.sh >/dev/null 2>&1 || true
    exit 0
  fi

  log "iteration $i / $MAX_ITERS (errors=$ERRORS, sorry warnings=$SORRIES, sorry tokens=$SORRY_LINES)"
  tag="$(printf '%03d' "$i")"
  snapshot="$(git stash create 2>/dev/null || true)"

  # The agent prompt: AGENT.md plus the current build status.
  prompt="$(cat cluster/AGENT.md)

---- CURRENT STATE (iteration $i) ----
lake build exit code: $BUILD_STATUS; errors: $ERRORS; 'sorry' warnings: $SORRIES; literal 'sorry' tokens: $SORRY_LINES
Axiom audit: $( [ "$AUDIT_OK" -eq 1 ] && echo clean || echo 'NOT clean (see cluster/last_audit.log)' )
Last build log (tail):
$(tail -n 60 cluster/last_build.log)
"
  # shellcheck disable=SC2086
  "$CLAUDE_BIN" -p "$prompt" $PERMISSION_MODE $CLAUDE_EXTRA \
    > "cluster/log/agent_${tag}.log" 2>&1
  agent_status=$?
  log "agent exit $agent_status (transcript cluster/log/agent_${tag}.log)"

  build_once
  audit_once
  cp cluster/last_build.log "cluster/log/build_${tag}.log"
  log "after iteration $i: build exit $BUILD_STATUS, errors=$ERRORS, sorry warnings=$SORRIES, sorry tokens=$SORRY_LINES, audit_ok=$AUDIT_OK"

  improved=0
  if [ "$ERRORS" -eq 0 ] && { [ "$SORRY_LINES" -lt "$best_sorries" ] || [ "$best_errors" -gt 0 ]; }; then
    improved=1
  fi

  if [ "$improved" -eq 1 ]; then
    best_sorries="$SORRY_LINES"; best_errors="$ERRORS"
    if [ "$NO_COMMIT" != "1" ]; then
      git add -A TwoTeacherCensus TwoTeacherCensus.lean PROOF_PLAN.md lakefile.lean 2>/dev/null
      git commit -q -m "cluster: iteration $i -- errors=$ERRORS, sorries=$SORRY_LINES, audit_ok=$AUDIT_OK" \
        -m "Automated proof pass; see cluster/log/agent_${tag}.log and cluster/log/build_${tag}.log." \
        2>/dev/null && log "committed iteration $i" || log "nothing to commit"
    fi
  elif [ "$ERRORS" -gt 0 ] && [ "$best_errors" -eq 0 ]; then
    # The agent broke the build; roll the tracked Lean files back to the last good state.
    log "build broken by iteration $i; reverting tracked files"
    git checkout -- TwoTeacherCensus TwoTeacherCensus.lean 2>/dev/null || true
    build_once
    audit_once
    log "after revert: errors=$ERRORS, sorry tokens=$SORRY_LINES"
  else
    log "no improvement in iteration $i"
  fi
done

log "stopped after $MAX_ITERS iterations: errors=$ERRORS, sorry tokens=$SORRY_LINES, audit_ok=$AUDIT_OK"
exit 1
