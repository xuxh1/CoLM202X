#!/bin/bash
# Make a V2 global continuation case in one call (does not submit):
#   create_clone of a finished global case -> bld/ source replaced by a
#   methane-paper-v2 commit (git archive, v2/ left out), keeping the source
#   case's define.h, Makeoptions link and every bld/run *.nml -> ifx cmake build
#   (log ends with MAKE_EXIT=<code>) -> copy_restart of <year0>-01-01 -> run
#   years -> DEF_METHANE overrides in bld/run/ch4_parameter.nml ->
#   CODE_COMMIT.txt. Only for changes that leave the carbon-pool initial state
#   alone (target doc section 1, item 5). Run on mgt01.
# Usage: [COLD=1] mk_global_cont.sh <source case> <dest case> <commit> <year0> <year1> "<note>" [KEY=VALUE ...]
#   case paths absolute or relative to cases/paper_v2; DEF_METHANE% keys go to
#   bld/run/ch4_parameter.nml, every other key to the &nl_colm group of the input
#   namelist. COLD=1 makes a cold-start case instead (no restart copied; submit
#   with lsf.sh submit <case> init so mkini runs first), for changes that do touch
#   the initial carbon pools.
set -euo pipefail
P=/share/home/dq076/mode/Methane/cases/paper_v2
abs() { case "$1" in /*) echo "$1" ;; *) echo "$P/$1" ;; esac; }
SRC=$(abs "${1:?source}"); DST=$(abs "${2:?dest}"); COMMIT=${3:?commit}; Y0=${4:?year0}; Y1=${5:?year1}
NOTE=${6:-}; shift 6
REPO=/share/home/dq076/mode/Methane/CoLM202X-paper-v2
FULL=$(git -C $REPO rev-parse "$COMMIT")
SN=$(basename "$SRC"); DN=$(basename "$DST")
[ -e "$DST" ] && { echo "ERROR: $DST exists" >&2; exit 1; }
[ "${COLD:-0}" = 1 ] || [ -d "$SRC/restart/$Y0-001-00000" ] || { echo "ERROR: no $SRC/restart/$Y0-001-00000" >&2; exit 1; }
mkdir -p "$(dirname "$DST")"

# 1. clone; the tool returns 0 even on failure, so check what it made
CL=$( cd $REPO/run/scripts && ./create_clone -s "$SRC" -d "$DST" 2>&1 ) || true
echo "$CL" > "$DST/clone.log" 2>/dev/null || { echo "ERROR: clone made no $DST" >&2; echo "$CL" >&2; exit 1; }
[ -f "$DST/input_$DN.nml" ] && [ -d "$DST/bld/main" ] || { echo "ERROR: clone incomplete, see $DST/clone.log" >&2; exit 1; }
[ -z "$(git -C $REPO status --short | grep -E '(^| )cases/' || true)" ] || { echo "ERROR: clone polluted $REPO" >&2; exit 1; }

# 2. source of the commit into bld/, source case's build settings and namelists back on top
[ -d "$DST/bld/build-dq076-ifx" ] && mv "$DST/bld/build-dq076-ifx" "$DST/bld/build-dq076-ifx.cloned"
git -C $REPO archive "$FULL" -- . ':(exclude)v2' | tar -x -C "$DST/bld"
cp -p "$SRC/bld/include/define.h" "$DST/bld/include/define.h"
ln -sfn Makeoptions.dq076-ifx "$DST/bld/include/Makeoptions"
( cd "$SRC/bld/run" && find . -name '*.nml' -print0 ) | while IFS= read -r -d '' f; do
  mkdir -p "$DST/bld/run/$(dirname "$f")"; cp -p "$SRC/bld/run/$f" "$DST/bld/run/$f"; done

# 3. build
LOG="$DST/build_$(date +%Y%m%d)_v2.log"
( set +u; source /share/home/dq076/software/intel-env-ifx && cd "$DST/bld" && ./cmake_build.sh ) > "$LOG" 2>&1 && rc=0 || rc=$?
echo "MAKE_EXIT=$rc" >> "$LOG"
[ "$rc" = 0 ] || { echo "ERROR: build failed, see $LOG" >&2; exit 1; }

# 4. restart of the start year, renamed to the new case by the tool (not for a cold start)
if [ "${COLD:-0}" != 1 ]; then
  ( cd $REPO/run/scripts && ./copy_restart -s "$SRC" -d "$DST" -o "$Y0-01-01" ) > "$DST/copy_restart.log" 2>&1 || true
  ls "$DST/restart/$Y0-001-00000/${DN}_restart_$Y0-001-00000"* > /dev/null 2>&1 || { echo "ERROR: no restart copied, see $DST/copy_restart.log" >&2; exit 1; }
fi

# 5. run years
sed -i -E "s/(DEF_simulation_time%start_year *= *)[0-9]{4}/\1$Y0/; s/(DEF_simulation_time%end_year *= *)[0-9]{4}/\1$Y1/" "$DST/input_$DN.nml"

# 6. DEF_METHANE overrides: replace in place, else add before the standalone '/'
#    closing &nl_colm_methane_parameter (the last DEF_METHANE% line can sit in
#    &nl_colm_tracer_forcing, which then rejects the key: V-26's first try)
CH4="$DST/bld/run/ch4_parameter.nml"; INP="$DST/input_$DN.nml"
for kv in "$@"; do
  k=${kv%%=*}; v=${kv#*=}; ke=$(printf '%s' "$k" | sed 's/[%]/\\%/g')
  case "$k" in DEF_METHANE%*) F=$CH4 ;; *) F=$INP ;; esac
  if grep -q "^ *$ke *=" "$F"; then sed -i "s|^\( *\)$ke *=.*|\1$k = $v|" "$F"
  elif [ "$F" = "$INP" ]; then
    n=$(grep -n '^ *&nl_colm *$' "$F" | head -1 | cut -d: -f1)
    [ -n "$n" ] || { echo "ERROR: no &nl_colm group in $F" >&2; exit 1; }
    sed -i "${n}a\\   $k = $v" "$F"
  else
    g=$(grep -n '^ *&nl_colm_methane_parameter' "$F" | head -1 | cut -d: -f1)
    n=$(awk -v g="$g" 'NR>g && /^ *\/ *$/ {print NR; exit}' "$F")
    [ -n "$g" ] && [ -n "$n" ] || { echo "ERROR: no &nl_colm_methane_parameter group in $F" >&2; exit 1; }
    sed -i "${n}i\\    $k = $v" "$F"
  fi
done

# 7. version record
cat > "$DST/CODE_COMMIT.txt" <<EOT
# 本 case 所用源码版本
repo    : $REPO
branch  : methane-paper-v2
commit  : $FULL
date    : $(git -C $REPO log -1 --format=%ad --date=iso "$FULL")
subject : $(git -C $REPO log -1 --format=%s "$FULL")
synced  : $(date '+%Y-%m-%d %H:%M') via create_clone（基线 ${SRC#/share/home/dq076/mode/Methane/cases/}）+ git archive 入 bld（v2/ 除外）
note    : define.h、Makeoptions 链与 bld/run 下全部 nml 取自基线；$([ "${COLD:-0}" = 1 ] && echo "冷启动（重做 mkini）跑" || echo "从基线 $Y0-001 重启续跑") $Y0–$Y1。$NOTE
EOT
echo "made $DST at $FULL, years $Y0-$Y1, overrides: $*"
diff "$SRC/bld/run/ch4_parameter.nml" "$CH4" || true
diff "$SRC/input_$SN.nml" "$DST/input_$DN.nml" || true
