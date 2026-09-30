#!/bin/bash
# Score V2 site tags as soon as they finish, so results never wait on a person.
# Reads v2/results/auto_tags.txt, one line per tag:
#   <case under cases/> <tag> main|val [eq-reference tag]
# When the tag's _conf/STATUS.txt says DONE and it has not been scored yet, runs
# on a compute node: summarize_tag.py (main) or summarize_val.py (val), then
# pathway_shares.py and gpp_check.py for main tags and eq_compare.py against the reference tag if
# one is given. Prints one line per scored tag (a Monitor event) and leaves a
# marker v2/results/scored/<tag>. Lines can be appended while it runs.
# Usage: auto_score.sh [node]   (loops until killed)
NODE=${1:-node112}
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
S=$V2/scripts
C=/share/home/dq076/mode/Methane/cases
mkdir -p $V2/results/scored
while true; do
  while read -r CASE TAG KIND REF; do
    [ -z "$TAG" ] && continue
    [ -f "$V2/results/scored/$TAG" ] && continue
    [ "$(head -1 $C/$CASE/sites_$TAG/_conf/STATUS.txt 2>/dev/null)" = DONE ] || continue
    if [ "$KIND" = val ]; then
      out=$($S/run_on_node.sh $NODE summarize_val.py "$CASE" "$TAG" 2>&1 | grep -m1 '^all 18')
    else
      out=$($S/run_on_node.sh $NODE summarize_tag.py "$CASE" "$TAG" 2>&1 | grep -m1 '^20 non-saline')
      $S/run_on_node.sh $NODE pathway_shares.py sites "$CASE" "$TAG" > $V2/results/pathway_shares_$TAG.log 2>&1
      $S/run_on_node.sh $NODE gpp_check.py "$CASE" "$TAG" > $V2/results/gpp_check_$TAG.log 2>&1
      out="$out | GPP ratio $(grep -o 'median ratio model/obs [0-9.]*' $V2/results/gpp_check_$TAG.log | awk '{print $4}')"
    fi
    eq=""
    if [ -n "$REF" ] && [ -f "$V2/results/scored/$REF" ]; then
      eq=$($S/run_on_node.sh $NODE eq_compare.py "$CASE" "$REF" "$TAG" 2>&1 | grep "^median |d|" | head -1)
    fi
    echo "$(date +%H:%M) scored $TAG: $out ${eq:+| vs $REF $eq}"
    date > "$V2/results/scored/$TAG"
  done < $V2/results/auto_tags.txt
  sleep 60
done
