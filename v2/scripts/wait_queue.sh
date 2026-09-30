#!/bin/bash
# Wait until N more site tags have been scored by the queue runner (lines
# appended to v2/results/queue_summary.tsv), then print the new lines' tag,
# subset and finish time. One call, so a background waiter needs no loop.
# Usage: wait_queue.sh <n>
N=${1:?n}
F=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2/results/queue_summary.tsv
N0=$(wc -l < "$F")
until [ "$(wc -l < "$F")" -ge $((N0 + N)) ]; do sleep 120; done
tail -n "$N" "$F" | cut -f1,4,5
