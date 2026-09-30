#!/bin/bash
# Wait for a V2 global continuation job, then run the standard evaluation suite
# on node114, in one call (so a background waiter needs no compound command).
# Prints the number of history files, the balance-error count in the log, and
# the result files written; if the last year is missing, the log tail instead.
# Usage: wait_eval_global.sh <case under cases/paper_v2> <jobid> <year0> <year1>
set -o pipefail
CASE=${1:?case}; JOB=${2:?jobid}; Y0=${3:?year0}; Y1=${4:?year1}
S=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/scripts
C=/share/home/dq076/mode/Methane/cases/paper_v2/$CASE
NAME=$(basename "$CASE")
"$S/wait_done.sh" global "$C" "$JOB"
echo "history files: $(ls "$C/history" | wc -l)"
echo "balance errors: $(grep -c 'balance error' "$C/log")"
if ls "$C/history/"*"_$Y1.nc" >/dev/null 2>&1; then
  "$S/global_suite.sh" "paper_v2/$(dirname "$CASE")" "$NAME" "$Y0" "$Y1" node114 >/dev/null 2>&1
  ls "$S/../results" | grep "$NAME"
else
  tail -5 "$C/log"
fi
