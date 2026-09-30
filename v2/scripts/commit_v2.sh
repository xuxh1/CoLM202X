#!/bin/bash
# Commit V2 records on methane-paper-v2 (never pushes). Paths are added with -f
# so the small result logs the ledger cites are kept although *.log is ignored.
# Usage: commit_v2.sh "<message>" <path relative to the worktree>...
set -eo pipefail
MSG=${1:?message}; shift
R=/share/home/dq076/mode/Methane/CoLM202X-paper-v2
[ "$(git -C $R rev-parse --abbrev-ref HEAD)" = methane-paper-v2 ] || { echo "not on methane-paper-v2" >&2; exit 1; }
git -C $R add -f -- "$@"
git -C $R commit -q -m "$MSG

git -C $R log --oneline -1
