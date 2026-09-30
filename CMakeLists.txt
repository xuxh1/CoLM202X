# CoLM202X -- CMake/Ninja front end for the Makefile.
#
# The Makefile stays the single source of truth. This file asks it (make -n -p)
# which compiler, flags and object lists include/Makeoptions and include/define.h
# resolve to, and hands the very same recipe to Ninja, which scans the Fortran
# module dependencies itself and therefore can build in parallel. No flag and no
# file name is duplicated here; add new sources to the Makefile as before.
#
#   ln -sfn Makeoptions.<stack> include/Makeoptions
#   cmake -S . -B build -G Ninja
#   cmake --build build -j4          # -> run/mksrfdata.x run/mkinidata.x run/colm.x
#   ./cmake_build.sh [-jN]           # the same, build dir chosen from the Makeoptions link
#
# Editing include/define.h or the Makefile re-runs the configure step on the
# next build. Use one build directory per compiler stack.
#
# Not covered (use make): CaMa_Flood builds and the `lib` target.

cmake_minimum_required(VERSION 3.20)

# ---------------------------------------------------------------------------
# 1. Ask the Makefile what it would do.
#    A one-target helper makefile echoes the expanded variables; `make -n`
#    prints those echo lines without running anything, `-p` adds the rule
#    database, from which the explicit `X.o: dir/Y.F90` overrides are read.
# ---------------------------------------------------------------------------
set(COLM_MAKE_VARS FF FOPTS LINK_FOPTS LDFLAGS NETCDF_INC VPATH CaMa
    OBJS_SHARED TRACER_INIT_CONFIG_OBJS TRACER_RUNTIME_CONFIG_OBJS
    TRACER_MKSRFDATA_CONFIG_OBJS OBJS_MKSRFDATA OBJS_BASIC OBJS_MKINIDATA
    OBJS_MAIN OBJS_POST1 OBJS_POST2 OBJS_POST3)

set(_colm_vars_mk "${CMAKE_BINARY_DIR}/colm_make_vars.mk")
set(_colm_mk "colm-cmake-vars:\n")
foreach(_v IN LISTS COLM_MAKE_VARS)
  string(APPEND _colm_mk "\t@echo COLMVAR ${_v} $(${_v})\n")
endforeach()
file(WRITE "${_colm_vars_mk}" "${_colm_mk}")

execute_process(
  COMMAND make -n -p -f Makefile -f "${_colm_vars_mk}" colm-cmake-vars
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE _colm_make_db
  ERROR_VARIABLE _colm_make_err
  RESULT_VARIABLE _colm_make_rc)
if(NOT _colm_make_rc EQUAL 0)
  message(FATAL_ERROR "make could not read the Makefile (is include/Makeoptions linked?):\n${_colm_make_err}")
endif()

set(_colm_make_db "\n${_colm_make_db}")   # so the first echo line also starts a line
foreach(_v IN LISTS COLM_MAKE_VARS)
  if(_colm_make_db MATCHES "\necho COLMVAR ${_v} ([^\n]*)")
    string(STRIP "${CMAKE_MATCH_1}" COLM_${_v})
  else()
    set(COLM_${_v} "")
  endif()
endforeach()

if(COLM_CaMa STREQUAL "YES")
  message(FATAL_ERROR "CaMa_Flood is on in include/define.h; CMakeLists.txt does not build CaMa, use make.")
endif()
if(NOT COLM_FF)
  message(FATAL_ERROR "FF is empty: include/Makeoptions did not define the compiler.")
endif()

set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS
  "${CMAKE_SOURCE_DIR}/Makefile" "${CMAKE_SOURCE_DIR}/include/Makeoptions"
  "${CMAKE_SOURCE_DIR}/include/define.h")

# ---------------------------------------------------------------------------
# 2. Compiler and flags: FF, FOPTS, LINK_FOPTS, LDFLAGS verbatim from make.
# ---------------------------------------------------------------------------
set(CMAKE_Fortran_COMPILER "${COLM_FF}" CACHE FILEPATH "FF from include/Makeoptions")
project(CoLM LANGUAGES Fortran)

separate_arguments(_colm_fopts      NATIVE_COMMAND "${COLM_FOPTS}")
separate_arguments(_colm_link_fopts NATIVE_COMMAND "${COLM_LINK_FOPTS}")
separate_arguments(_colm_ldflags    NATIVE_COMMAND "${COLM_LDFLAGS}")

# Keep the build-type flag sets empty so nothing is added on top of FOPTS.
set(CMAKE_Fortran_FLAGS "")
foreach(_c RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
  set(CMAKE_Fortran_FLAGS_${_c} "")
endforeach()
add_compile_options(${_colm_fopts})

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/modules")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/run")
include_directories("${CMAKE_SOURCE_DIR}/include" "${CMAKE_BINARY_DIR}/modules")
if(COLM_NETCDF_INC)
  include_directories("${COLM_NETCDF_INC}")
endif()

# ---------------------------------------------------------------------------
# 3. Object name -> source file, the way make's VPATH would resolve it,
#    plus the explicit `X.o: dir/Y.F90` rules (extends/interception).
# ---------------------------------------------------------------------------
string(REPLACE ":" ";" _colm_vpath "${COLM_VPATH}")
foreach(_d IN LISTS _colm_vpath)
  string(STRIP "${_d}" _d)
  if(_d STREQUAL "" OR _d STREQUAL ".bld")
    continue()
  endif()
  file(GLOB _f90 RELATIVE "${CMAKE_SOURCE_DIR}" CONFIGURE_DEPENDS "${CMAKE_SOURCE_DIR}/${_d}/*.F90")
  foreach(_f IN LISTS _f90)
    get_filename_component(_n "${_f}" NAME_WE)
    if(NOT DEFINED COLM_SRC_${_n})        # first VPATH hit wins, as in make
      set(COLM_SRC_${_n} "${_f}")
    endif()
  endforeach()
endforeach()

string(REGEX MATCHALL "\n[A-Za-z0-9_]+\\.o: [^ \n]*/[^ \n]*\\.F90" _colm_overrides "${_colm_make_db}")
foreach(_r IN LISTS _colm_overrides)
  string(REGEX REPLACE "\n([A-Za-z0-9_]+)\\.o: (.*)" "\\1" _n "${_r}")
  string(REGEX REPLACE "\n([A-Za-z0-9_]+)\\.o: (.*)" "\\2" _f "${_r}")
  set(COLM_OVERRIDE_${_n} "${_f}")
endforeach()

function(colm_sources out_var)
  set(_srcs)
  foreach(_o IN LISTS ARGN)
    string(REGEX REPLACE "\\.o$" "" _n "${_o}")
    if(DEFINED COLM_OVERRIDE_${_n})
      list(APPEND _srcs "${COLM_OVERRIDE_${_n}}")
    elseif(DEFINED COLM_SRC_${_n})
      list(APPEND _srcs "${COLM_SRC_${_n}}")
    else()
      message(FATAL_ERROR "no .F90 found for ${_o} in VPATH (${COLM_VPATH})")
    endif()
  endforeach()
  set(${out_var} "${_srcs}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# 4. One object library per Makefile stage; a source that appears in several
#    lists is compiled once, in the earliest stage. Stage order follows the
#    Makefile's link lines, so module files exist before they are used.
#    Executables take the objects explicitly ($<TARGET_OBJECTS>), so what the
#    stage-to-stage links propagate is module information only.
# ---------------------------------------------------------------------------
set(_colm_seen "")
macro(colm_stage name)
  set(_objs ${ARGN})
  if(_colm_seen)
    list(REMOVE_ITEM _objs ${_colm_seen})
  endif()
  if(_objs)
    list(REMOVE_DUPLICATES _objs)
    colm_sources(_srcs ${_objs})
    add_library(${name} OBJECT ${_srcs})
    list(APPEND _colm_seen ${_objs})
  endif()
endmacro()

# Linking one object library to another carries no objects, only the Fortran
# module provenance Ninja needs to order the compiles across stages.
function(colm_after target)
  if(TARGET ${target})
    foreach(_d IN LISTS ARGN)
      if(TARGET ${_d})
        target_link_libraries(${target} PUBLIC ${_d})
      endif()
    endforeach()
  endif()
endfunction()

function(colm_executable name)
  set(_objs)
  foreach(_l IN LISTS ARGN)
    if(TARGET ${_l})
      list(APPEND _objs "$<TARGET_OBJECTS:${_l}>")
    endif()
  endforeach()
  add_executable(${name} ${_objs})
  target_link_options(${name} PRIVATE ${_colm_link_fopts})
  target_link_libraries(${name} PRIVATE ${_colm_ldflags})
endfunction()

foreach(_v OBJS_SHARED TRACER_INIT_CONFIG_OBJS TRACER_RUNTIME_CONFIG_OBJS
           TRACER_MKSRFDATA_CONFIG_OBJS OBJS_MKSRFDATA OBJS_BASIC OBJS_MKINIDATA
           OBJS_MAIN OBJS_POST1 OBJS_POST2 OBJS_POST3)
  separate_arguments(_colm_${_v} NATIVE_COMMAND "${COLM_${_v}}")
endforeach()

colm_stage(colm_shared        ${_colm_OBJS_SHARED})
colm_stage(colm_tracer_config ${_colm_TRACER_INIT_CONFIG_OBJS}
                              ${_colm_TRACER_RUNTIME_CONFIG_OBJS}
                              ${_colm_TRACER_MKSRFDATA_CONFIG_OBJS})
colm_stage(colm_mksrfdata     ${_colm_OBJS_MKSRFDATA})
colm_stage(colm_basic         ${_colm_OBJS_BASIC})
colm_stage(colm_mkinidata     ${_colm_OBJS_MKINIDATA})
colm_stage(colm_main          ${_colm_OBJS_MAIN})
colm_stage(colm_post1         ${_colm_OBJS_POST1})
colm_stage(colm_post2         ${_colm_OBJS_POST2})
colm_stage(colm_post3         ${_colm_OBJS_POST3})

colm_after(colm_tracer_config colm_shared)
colm_after(colm_mksrfdata     colm_shared colm_tracer_config)
colm_after(colm_basic         colm_shared colm_tracer_config)
colm_after(colm_mkinidata     colm_shared colm_tracer_config colm_basic)
colm_after(colm_main          colm_shared colm_tracer_config colm_basic)
colm_after(colm_post1         colm_shared)
colm_after(colm_post2         colm_shared)
colm_after(colm_post3         colm_shared)

# The three programs the project runs. The tracer config objects go into all
# three (make links a subset into each); they use only shared and tracer
# config modules, so the extra objects are unreferenced and change nothing.
colm_executable(mksrfdata.x colm_shared colm_tracer_config colm_mksrfdata)
colm_executable(mkinidata.x colm_shared colm_tracer_config colm_basic colm_mkinidata)
colm_executable(colm.x      colm_shared colm_tracer_config colm_basic colm_main)

# Post-processing tools, on demand: cmake --build build --target hist_concatenate.x
foreach(_p "hist_concatenate.x;colm_post1" "post_vector2grid.x;colm_post2" "srfdata_concatenate.x;colm_post3")
  list(GET _p 0 _exe)
  list(GET _p 1 _lib)
  if(TARGET ${_lib})
    colm_executable(${_exe} colm_shared ${_lib})
    set_target_properties(${_exe} ${_lib} PROPERTIES EXCLUDE_FROM_ALL ON)
  endif()
endforeach()

message(STATUS "CoLM: FF      = ${COLM_FF}")
message(STATUS "CoLM: FOPTS   = ${COLM_FOPTS}")
message(STATUS "CoLM: LDFLAGS = ${COLM_LDFLAGS}")
list(LENGTH _colm_seen _n_objs)
message(STATUS "CoLM: ${_n_objs} objects from the Makefile lists")
