#!/bin/bash
# Launch one V2 site tag in one call: write its namelists from a V2 config,
# submit it to LSF, and register it for auto_score.sh.
# Usage: launch.sh <case under cases/> <config under v2/config> <tag> main|val [spinup] [eq-reference tag]
#   main uses the config's LIST_sites.csv, SITE_PARAMS.csv (and SITE_MAIN.csv); val uses the
#   validation list v2/config/LIST_sites_val23.csv and SITE_PARAMS_val23_c13.csv.
set -eo pipefail
CASE=${1:?case}; CFGN=${2:?config}; TAG=${3:?tag}; KIND=${4:?main|val}; SPIN=${5:-30}; REF=${6:-}
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
S=$V2/scripts; CF=$V2/config
if [ "$KIND" = val ]; then
  LIST=$CF/LIST_sites_val23.csv; SP=$CF/SITE_PARAMS_val23_c13.csv; SM=$CF/SITE_MAIN_val23.csv
  # a config may carry its own validation site parameters (V-47)
  [ -f "$CF/$CFGN/SITE_PARAMS_val23.csv" ] && SP=$CF/$CFGN/SITE_PARAMS_val23.csv
else
  LIST=$CF/$CFGN/LIST_sites.csv; SP=$CF/$CFGN/SITE_PARAMS.csv; SM=$CF/$CFGN/SITE_MAIN.csv
fi
# per-site main-namelist keys (C-17c peat share) only for configs that carry SITE_MAIN.csv,
# so older configs and trees never see a key their namelist lacks
if [ -f "$CF/$CFGN/SITE_MAIN.csv" ] && [ -f "$SM" ]; then :; else SM=""; fi
$S/prep_tag.sh "$CASE" "$TAG" "$CFGN" "$LIST" "$SPIN" "$SP" "$SM" | tail -1 | sed 's#.*/cases/#prepared cases/#'
$S/lsf.sh sites "$CASE" "$TAG"
echo "$CASE $TAG $KIND $REF" >> $V2/results/auto_tags.txt
