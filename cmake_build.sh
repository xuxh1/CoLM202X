#!/bin/bash
# CoLM CMake/Ninja build: configure once per stack, then build.
#   ./cmake_build.sh                       # -j4 on mgt01, idle cores elsewhere
#   ./cmake_build.sh -j40
#   ./cmake_build.sh --target hist_concatenate.x
# The build directory is build-<stack>, <stack> read from the include/Makeoptions
# link (Makeoptions.dq089-gnu -> build-dq089-gnu), so relinking Makeoptions
# switches directories and never mixes two stacks in one cache. Source the
# matching environment (gnu-env or intel-env-ifx) before running, as for make.
set -e
cd "$(dirname "$0")"
CMAKE=/share/home/dq076/software/miniconda3/envs/py311/bin/cmake
NINJA=/share/home/dq076/software/miniconda3/envs/py311/bin/ninja
stack=$(readlink include/Makeoptions | sed 's/^Makeoptions\.//')
if [ -z "$stack" ]; then
    echo "include/Makeoptions is not a symlink; run: ln -sfn Makeoptions.<stack> include/Makeoptions" >&2
    exit 1
fi
bld=build-$stack
[ -f "$bld/build.ninja" ] || "$CMAKE" -S . -B "$bld" -G Ninja -DCMAKE_MAKE_PROGRAM="$NINJA"
case " $* " in
    *" -j"*) ;;
    *)  if [ "$(hostname)" = mgt01 ]; then
            j=4                                   # per-user cgroup: 4 cores / 8 GB
        else
            load=$(cut -d. -f1 /proc/loadavg)     # 1-minute load, integer part
            j=$(( $(nproc) - load )); [ "$j" -ge 1 ] || j=1
        fi
        set -- -j"$j" "$@" ;;
esac
exec "$CMAKE" --build "$bld" "$@"
