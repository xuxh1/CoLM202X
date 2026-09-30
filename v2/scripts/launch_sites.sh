#!/bin/bash
# Start a site tag on a compute node through the project driver (no copies of
# scripts/sites in the case). Usage: launch_sites.sh <node> <case under cases/> <tag> [workers] [run_sites args...]
set -eo pipefail
NODE=${1:?node}; CASE=${2:?case}; TAG=${3:?tag}; W=${4:-24}; shift 4 2>/dev/null || shift $#
cd /share/home/dq076/mode/Methane
WORKERS=$W scripts/sites/run_sites_direct.sh "$NODE" "$CASE" "$TAG" "$@" 2>&1 | grep -E "started|ERROR|refus"
