#!/bin/bash
# Build a single-point tree for paper V2 from a commit of methane-paper-v2.
# Usage: build_site_tree.sh <commit> <case dir under cases/paper_v2, or absolute> "<note>"
# Steps: git archive (without v2/) -> define.h to single point -> Makeoptions
# link to the ifx stack -> cmake build (log ends with MAKE_EXIT=<code>) ->
# CODE_COMMIT.txt -> config snapshot. Run on mgt01 (links MKL LAPACK).
set -euo pipefail
COMMIT=${1:?commit}; CASE=${2:?case dir}; NOTE=${3:-}
case "$CASE" in /*) ;; *) CASE=/share/home/dq076/mode/Methane/cases/paper_v2/$CASE ;; esac   # relative to cases/paper_v2
REPO=/share/home/dq076/mode/Methane/CoLM202X-paper-v2
FULL=$(git -C $REPO rev-parse "$COMMIT")
[ -e "$CASE" ] && { echo "ERROR: $CASE exists" >&2; exit 1; }
mkdir -p "$CASE"
git -C $REPO archive "$FULL" -- . ':(exclude)v2' | tar -x -C "$CASE"
sed -i 's/^#define GRIDBASED/#undef GRIDBASED/; s/^#undef SinglePoint/#define SinglePoint/' "$CASE/include/define.h"
ln -sfn Makeoptions.dq076-ifx "$CASE/include/Makeoptions"
LOG="$CASE/build_$(date +%Y%m%d)_v2.log"
( set +u; source /share/home/dq076/software/intel-env-ifx && cd "$CASE" && ./cmake_build.sh ) > "$LOG" 2>&1 && rc=0 || rc=$?
echo "MAKE_EXIT=$rc" >> "$LOG"
cat > "$CASE/CODE_COMMIT.txt" <<EOT
# 本 case 所用源码版本
repo    : $REPO
branch  : methane-paper-v2
commit  : $FULL
date    : $(git -C $REPO log -1 --format=%ad --date=iso "$FULL")
subject : $(git -C $REPO log -1 --format=%s "$FULL")
synced  : $(date '+%Y-%m-%d %H:%M') via git archive (v2/ excluded)
note    : define.h 由网格改单点（GRIDBASED 关、SinglePoint 开），其余逐字同源码树。$NOTE
EOT
/share/home/dq076/mode/Methane/scripts/snapshot_case_config.sh "$CASE" > /dev/null
echo "built $CASE at $FULL, MAKE_EXIT=$rc"
