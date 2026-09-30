#!/bin/bash
# Bring an existing V2 global case to a newer methane-paper-v2 commit before
# its model run is submitted (for example after its mksrf ran from an older
# build): the source of the commit goes into bld/ over the old one (v2/ left
# out; define.h, the Makeoptions link and every bld/run *.nml kept), the ifx
# cmake build runs again (log ends with MAKE_EXIT=<code>), CODE_COMMIT.txt is
# rewritten, and KEY=VALUE overrides go into the namelists as in
# mk_global_cont.sh. Refuses while an LSF job of the case is pending or
# running. Run on mgt01.
# Usage: resync_case.sh <case dir> <commit> "<note>" [KEY=VALUE ...]
set -euo pipefail
P=/share/home/dq076/mode/Methane/cases/paper_v2
REPO=/share/home/dq076/mode/Methane/CoLM202X-paper-v2
C=$(case "${1:?case}" in /*) echo "$1" ;; *) echo "$P/$1" ;; esac); COMMIT=${2:?commit}; NOTE=${3:-}; shift 3
N=$(basename "$C"); FULL=$(git -C $REPO rev-parse "$COMMIT")
set +eu; source /share/lsf-suite/lsf/.conf_tmpl/profile.lsf > /dev/null 2>&1; set -eu
if bjobs -noheader -o "stat job_name" 2>/dev/null | grep -E '^(PEND|RUN) ' | awk '{print $2}' | grep -qx "$N"; then
  echo "ERROR: an LSF job named $N is pending or running" >&2; exit 1
fi
TMP=$(mktemp -d "$C/.resync.XXXX")
cp -p "$C/bld/include/define.h" "$TMP/define.h"
( cd "$C/bld/run" && find . -name '*.nml' -print0 ) | while IFS= read -r -d '' f; do
  mkdir -p "$TMP/run/$(dirname "$f")"; cp -p "$C/bld/run/$f" "$TMP/run/$f"; done
git -C $REPO archive "$FULL" -- . ':(exclude)v2' | tar -x -C "$C/bld"
cp -p "$TMP/define.h" "$C/bld/include/define.h"
ln -sfn Makeoptions.dq076-ifx "$C/bld/include/Makeoptions"
( cd "$TMP/run" && find . -name '*.nml' -print0 ) | while IFS= read -r -d '' f; do cp -p "$TMP/run/$f" "$C/bld/run/$f"; done
LOG="$C/build_$(date +%Y%m%d_%H%M)_resync.log"
( set +u; source /share/home/dq076/software/intel-env-ifx && cd "$C/bld" && ./cmake_build.sh ) > "$LOG" 2>&1 && rc=0 || rc=$?
echo "MAKE_EXIT=$rc" >> "$LOG"
[ "$rc" = 0 ] || { echo "ERROR: build failed, see $LOG" >&2; exit 1; }
CH4="$C/bld/run/ch4_parameter.nml"; INP="$C/input_$N.nml"
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}; ke=$(printf '%s' "$k" | sed 's/[%]/\\%/g')
  case "$k" in DEF_METHANE%*) F=$CH4 ;; *) F=$INP ;; esac
  if grep -q "^ *$ke *=" "$F"; then sed -i "s|^\( *\)$ke *=.*|\1$k = $v|" "$F"
  elif [ "$F" = "$INP" ]; then
    n=$(grep -n '^ *&nl_colm *$' "$F" | head -1 | cut -d: -f1); sed -i "${n}a\\   $k = $v" "$F"
  else
    g=$(grep -n '^ *&nl_colm_methane_parameter' "$F" | head -1 | cut -d: -f1)
    n=$(awk -v g="$g" 'NR>g && /^ *\/ *$/ {print NR; exit}' "$F"); sed -i "${n}i\\    $k = $v" "$F"
  fi
done
OLD=$(grep '^commit' "$C/CODE_COMMIT.txt" | awk '{print $3}')
cat > "$C/CODE_COMMIT.txt" <<EOT
# 本 case 所用源码版本
repo    : $REPO
branch  : methane-paper-v2
commit  : $FULL
date    : $(git -C $REPO log -1 --format=%ad --date=iso "$FULL")
subject : $(git -C $REPO log -1 --format=%s "$FULL")
synced  : $(date '+%Y-%m-%d %H:%M') via resync_case.sh（原 ${OLD:0:8}，bld 源码换为本 commit，define.h、Makeoptions 链与 bld/run 下全部 nml 保留）
note    : $NOTE
EOT
rm -rf "$TMP"
echo "resynced $C to ${FULL:0:8}; overrides: $*"
