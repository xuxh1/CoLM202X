#!/bin/bash
# LSF from mgt01 for the V2 global cases.
# Usage: lsf.sh status <jobid...> | lsf.sh top <jobid> | lsf.sh submit <case dir> [init] [after <jobid>] | lsf.sh kill <jobid>
#        lsf.sh sites <case under cases/> <tag> [run_sites args...]
# submit: pre-creates restart/<year>-001-00000 for every simulated year and the one after
# (ifx mkdir can fail, V-3), snapshots the configuration, checks the generation, then
# bsub's case.submit (after init.submit when 'init' is given). 'after <jobid>' holds the
# first job until that job has ended (needed only to keep two mksrf jobs apart, iron rule 8).
set -eo pipefail
set +e; source /share/lsf-suite/lsf/.conf_tmpl/profile.lsf; set -e   # the profile returns non-zero off Cray hosts
R=/share/home/dq076/mode/Methane
case "${1:?status|submit|kill}" in
  status) shift; bjobs -a -o "jobid job_name stat start_time finish_time" -noheader "$@" ;;
  kill)   bkill "${2:?jobid}" ;;
  top)    btop "${2:?jobid}" ;;   # move one of my pending jobs to the front of my queue
  sites)
    # A site tag on one LSF node through scripts/sites/submit_sites.lsf, with
    # CASE, TAG and EXTRA filled in; the filled copy is kept in the tag's _conf/.
    CASE=${2:?case}; TAG=${3:?tag}; shift 3; EXTRA="$*"
    D=$R/cases/$CASE/sites_$TAG/_conf; [ -d "$D" ] || { echo "no $D; run prep_tag.sh first" >&2; exit 1; }
    mkdir -p $R/outputs/logs/lsf
    sed -e "s#^CASE=.*#CASE=$CASE#" -e "s#^TAG=.*#TAG=$TAG#" -e "s#^EXTRA=.*#EXTRA=\"$EXTRA\"#" \
        -e "s#^\#BSUB -J .*#\#BSUB -J v2_$TAG#" $R/scripts/sites/submit_sites.lsf > $D/submit_sites.lsf
    bsub < $D/submit_sites.lsf ;;
  submit)
    C=$(readlink -f "${2:?case dir}"); N=$(basename "$C"); shift 2; INIT=""; W=()
    while [ $# -gt 0 ]; do case "$1" in init) INIT=1 ;; after) W=(-w "ended($2)"); shift ;; esac; shift; done
    y0=$(grep -m1 'start_year' "$C/input_$N.nml" | grep -o '[0-9]\{4\}'); y1=$(grep -m1 'end_year' "$C/input_$N.nml" | grep -o '[0-9]\{4\}')
    mkdir -p "$C/restart/const"; for y in $(seq $y0 $((y1+1))); do mkdir -p "$C/restart/$y-001-00000"; done
    ( cd $R && scripts/snapshot_case_config.sh "${C#$R/}" > /dev/null 2>&1 ) || echo "snapshot failed" >&2
    ( cd $R && scripts/check_case_versions.sh "$(basename $(dirname $C))" | tail -1 ) || true
    cd "$C"
    if [ -n "$INIT" ]; then
      J=$(bsub "${W[@]}" < init.submit | grep -o '<[0-9]*>' | tr -d '<>'); echo "init $J"
      sed "s/^#BSUB -J \(.*\)$/#BSUB -J \1\n#BSUB -w done($J)/" case.submit > case_dep.submit; bsub < case_dep.submit
    else
      bsub "${W[@]}" < case.submit
    fi ;;
esac
