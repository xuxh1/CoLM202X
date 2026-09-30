#!/bin/bash
# Block until a V2 run has ended, then print one line (for a background waiter).
# Usage: wait_done.sh site <case under cases/> <tag>     ends when _conf/STATUS.txt leaves RUNNING
#        wait_done.sh global <case dir> <jobid>          ends when the LSF job leaves PEND/RUN
set -u
set +e; source /share/lsf-suite/lsf/.conf_tmpl/profile.lsf > /dev/null 2>&1; set -e
R=/share/home/dq076/mode/Methane/cases
case "${1:?site|global}" in
  site)
    S=$R/$2/sites_$3/_conf/STATUS.txt
    while [ "$(head -1 $S 2>/dev/null)" = RUNNING ] || [ ! -f "$S" ]; do sleep 300; done
    echo "$3: $(head -1 $S) $(sed -n 4p $S)" ;;
  global)
    while bjobs -noheader -o stat "$3" 2>/dev/null | grep -qE 'PEND|RUN'; do sleep 300; done
    echo "$(basename $2): job $3 $(bjobs -a -noheader -o stat $3 2>/dev/null); history: $(ls $2/history 2>/dev/null | grep -o '_20[0-9][0-9]' | sort -u | tr -d '_' | tr '\n' ' ')" ;;
esac
