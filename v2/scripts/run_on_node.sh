#!/bin/bash
# Run a V2 analysis script on a compute node (never on mgt01), with HDF5 file
# locking off so a history file the model has open is never locked (v2budget.py).
# Usage: run_on_node.sh <node> <script.py under v2/scripts> [args...]
#   V2_YEARS and V2_MONTHS in the environment are passed through (model_per_area.py, area_bands.py).
set -eo pipefail
NODE=${1:?node}; SCRIPT=${2:?script}; shift 2
V2=/share/home/dq076/mode/Methane/CoLM202X-paper-v2/v2
PY=/share/home/dq076/software/miniconda3/envs/py311/bin/python
ARGS=$(printf ' %q' "$@")
ssh -o ConnectTimeout=5 -o BatchMode=yes "$NODE" \
  "export HDF5_USE_FILE_LOCKING=FALSE V2_YEARS='${V2_YEARS:-}' V2_MONTHS='${V2_MONTHS:-}'; [ -z \"\$V2_YEARS\" ] && unset V2_YEARS; [ -z \"\$V2_MONTHS\" ] && unset V2_MONTHS; cd $V2 && $PY -W ignore scripts/$SCRIPT$ARGS" 2>&1 | grep -v getfattr
