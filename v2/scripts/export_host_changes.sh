#!/bin/bash
# Export the host-model (CoLM outside main/TRACER) part of every V2 code commit
# as its own patch, so host changes can be kept, applied or dropped apart from
# the methane module (user 2026-09-28: host code is not the user's to change or
# submit). Writes v2/host_changes/patches/NN_<hash>.patch (host paths only),
# v2/host_changes/host_cumulative.patch (all host changes since the V2 base) and
# v2/host_changes/commits.tsv (hash, date, host files, methane files in the same
# commit, subject). The hand-kept register is v2/host_changes/register.tsv.
# Usage: export_host_changes.sh
set -euo pipefail
R=/share/home/dq076/mode/Methane/CoLM202X-paper-v2
O=$R/v2/host_changes
BASE=$(git -C $R merge-base HEAD methane-paper)
HOST=(-- . ':(exclude)v2' ':(exclude)main/TRACER')
mkdir -p $O/patches
rm -f $O/patches/*.patch
printf 'n\tcommit\tdate\thost_files\tmethane_files_same_commit\tsubject\n' > $O/commits.tsv
n=0
for h in $(git -C $R log --reverse --format=%h $BASE..HEAD "${HOST[@]}"); do
  n=$((n+1)); nn=$(printf '%02d' $n)
  git -C $R format-patch -1 --stdout $h "${HOST[@]}" > $O/patches/${nn}_$h.patch
  hf=$(git -C $R show --name-only --format= $h "${HOST[@]}" | tr '\n' ' ' | sed 's/ $//')
  tf=$(git -C $R show --name-only --format= $h -- main/TRACER | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' $nn $h "$(git -C $R log -1 --format=%ad --date=short $h)" "$hf" "$tf" \
    "$(git -C $R log -1 --format=%s $h)" >> $O/commits.tsv
done
git -C $R diff $BASE HEAD "${HOST[@]}" > $O/host_cumulative.patch
echo "base $BASE; $n host commits; cumulative patch $(wc -l < $O/host_cumulative.patch) lines"
