#!/bin/bash
# One-screen status of V2 work: every tag in v2/results/auto_tags.txt (sites
# complete, STATUS, scored or not) and the given LSF jobs.
# Usage: status.sh [jobid ...]
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
C=/share/home/dq076/mode/Methane/cases
while read -r CASE TAG KIND REF; do
  [ -z "$TAG" ] && continue
  n=$(ls $C/$CASE/sites_$TAG/*/.complete 2>/dev/null | wc -l)
  st=$(head -1 $C/$CASE/sites_$TAG/_conf/STATUS.txt 2>/dev/null)
  sc=$([ -f $V2/results/scored/$TAG ] && echo scored || echo -)
  printf '%-22s %-5s %2s/23 %-8s %s\n' "$TAG" "$KIND" "$n" "${st:-?}" "$sc"
done < $V2/results/auto_tags.txt
if [ $# -gt 0 ]; then
  set +e; source /share/lsf-suite/lsf/.conf_tmpl/profile.lsf > /dev/null 2>&1; set -e
  bjobs -a -o "jobid job_name stat start_time finish_time" -noheader "$@"
fi
date +%H:%M
