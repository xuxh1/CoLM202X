#!/bin/bash
# All global checks of one V2 case in one call, each log into v2/results/:
# global_eval (totals, bands, area), region_eval (18 GCP regions), rice_regions,
# lake_regions, flood_regions, country_eval. Runs on a compute node.
# Usage: global_suite.sh <version dir under cases/> <case name> <y0> <y1> [node]
set -eo pipefail
VER=${1:?version}; N=${2:?case}; Y0=${3:?y0}; Y1=${4:?y1}; NODE=${5:-node111}
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
for s in global_eval region_eval rice_regions lake_regions flood_regions country_eval; do
  $V2/scripts/run_on_node.sh $NODE $s.py "$VER" "$N" "$Y0" "$Y1" > $V2/results/${s}_${N}_$Y0-$Y1.log 2>&1 || echo "$s failed"
done
grep -h -E '^E_|^# regions|^# median|^sum of' $V2/results/global_eval_${N}_$Y0-$Y1.log $V2/results/region_eval_${N}_$Y0-$Y1.log
grep -h -E "^# inside range" $V2/results/country_eval_${N}_$Y0-$Y1.log
