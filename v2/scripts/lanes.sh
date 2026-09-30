#!/bin/bash
# Two-lane check of the V2 loop. User 2026-09-30: three lanes spend tokens too
# fast; run two lanes, site (with literature folded into it) and global, and wait
# for the site results instead of filling every slot. Prints one line per lane
# and exits 1 only when a lane has nothing at all, so the Stop hook refuses to end
# a turn only when a lane has gone completely idle.
#   site:   pending + running rows of v2/queue.tsv, or rows of v2/三线.tsv with
#           status 进行中 in the 站点, 文献 or 实现 lanes; idle when none
#   global: running + pending g2_* LSF jobs of dq076; idle when none
# The three-lane version (site >= 4, global >= 3, literature and implementation
# each busy, backlog >= 2 per lane) is in git history before this change.
# Usage: lanes.sh
V=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
idle=0
s=$(awk -F'\t' 'NR>1 && ($2=="pending" || $2=="running")' $V/queue.tsv | wc -l)
w=$(awk -F'\t' '($1=="站点" || $1=="文献" || $1=="实现") && $3=="进行中"' $V/三线.tsv | wc -l)
if [ "$s" -lt 1 ] && [ "$w" -lt 1 ]; then
  echo "site+literature: IDLE (no site tag and no study in progress)"; idle=1
else
  echo "site+literature: ok ($s site tags, $w studies in progress)"
fi
source /share/lsf-suite/lsf/.conf_tmpl/profile.lsf > /dev/null 2>&1 || true
g=$(bjobs -w -noheader 2>/dev/null | awk '($3=="RUN" || $3=="PEND") && $7 ~ /^g2_/' | wc -l)
[ "$g" -lt 1 ] && { echo "global: IDLE (no g2 job)"; idle=1; } || echo "global: ok ($g running or pending)"
exit $idle
