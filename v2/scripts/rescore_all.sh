#!/bin/bash
# Rescore every V2 site tag that has a summary in results/pre_d17/ (case taken
# from the summary header) with the current summary scripts (D-17 gate), and
# print old -> new main lines. Run on a compute node.
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
PY=/share/home/dq076/software/miniconda3/envs/py311/bin/python
export HDF5_USE_FILE_LOCKING=FALSE
for f in $V2/results/pre_d17/*_summary.txt; do
  tag=$(basename "$f" _summary.txt)
  case "$tag" in *rice*|*twt*) continue;; esac
  case=$(head -1 "$f" | sed -n 's/^# [^ ]* (\([^)]*\)).*/\1/p')
  [ -n "$case" ] || { echo "$tag no case in header"; continue; }
  if head -1 "$f" | grep -q "validation"; then s=summarize_val.py; k='^all '; else s=summarize_tag.py; k='^20 non-saline'; fi
  $PY $V2/scripts/$s "$case" "$tag" > /dev/null 2>&1 || { echo "$tag FAILED"; continue; }
  o=$(grep -m1 "$k" "$f" | sed 's/.*: *//')
  n=$(grep -m1 "$k" $V2/results/${tag}_summary.txt | sed 's/.*: *//')
  printf '%-18s old %-28s new %s\n' "$tag" "$o" "$n"
done
