###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2018 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find PAMPA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PAMPA
#               [REQUIRED]             # Fail with error if pampa is not found
#              )
#
#  PAMPA depends on the following libraries:
#   - MPI
#   - PTSCOTCH
#
# This module finds headers and pampa library.
# Results are reported in variables:
#  PAMPA_FOUND             - True if headers and requested libraries were found
#  PAMPA_CFLAGS_OTHER      - pampa compiler flags without headers paths
#  PAMPA_LDFLAGS_OTHER     - pampa linker flags without libraries
#  PAMPA_INCLUDE_DIRS      - pampa include directories
#  PAMPA_LIBRARY_DIRS      - pampa link directories
#  PAMPA_LIBRARIES         - pampa libraries to be linked (absolute path)
#  PAMPA_CFLAGS_OTHER_DEP  - pampa + dependencies compiler flags without headers paths
#  PAMPA_LDFLAGS_OTHER_DEP - pampa + dependencies linker flags without libraries
#  PAMPA_INCLUDE_DIRS_DEP  - pampa + dependencies include directories
#  PAMPA_LIBRARY_DIRS_DEP  - pampa + dependencies link directories
#  PAMPA_LIBRARIES_DEP     - pampa + dependencies libraries
#  PAMPA_INTSIZE         - Number of octets occupied by a PAMPA_Num
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPAMPA=path/to/pampa):
#  PAMPA_DIR             - Where to find the base directory of pampa
#  PAMPA_INCDIR          - Where to find the header files
#  PAMPA_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PAMPA_DIR, PAMPA_INCDIR, PAMPA_LIBDIR

#=============================================================================
# Copyright 2012-2018 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2018 Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)

# Common macros to use in finds
include(FindInit)

if (NOT PAMPA_FOUND)
  set(PAMPA_DIR "" CACHE PATH "Installation directory of PAMPA library")
  if (NOT PAMPA_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PAMPA_DIR, has been set to specify the install directory of PAMPA")
  endif()
endif()

# PAMPA depends on MPI, try to find it
if (PAMPA_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# PAMPA depends on PAMPA, try to find it
if (PAMPA_FIND_REQUIRED)
  find_package(PTSCOTCH REQUIRED)
else()
  find_package(PTSCOTCH)
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_PAMPA_DIR "$ENV{PAMPA_DIR}")
set(ENV_PAMPA_INCDIR "$ENV{PAMPA_INCDIR}")
if(ENV_PAMPA_INCDIR)
  list(APPEND _inc_env "${ENV_PAMPA_INCDIR}")
elseif(ENV_PAMPA_DIR)
  list(APPEND _inc_env "${ENV_PAMPA_DIR}")
  list(APPEND _inc_env "${ENV_PAMPA_DIR}/include")
  list(APPEND _inc_env "${ENV_PAMPA_DIR}/include/ptscotch")
else()
  if(WIN32)
    string(REPLACE ":" ";" _inc_env "$ENV{INCLUDE}")
  else()
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{C_INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{CPATH}")
    list(APPEND _inc_env "${_path_env}")
    string(REPLACE ":" ";" _path_env "$ENV{INCLUDE_PATH}")
    list(APPEND _inc_env "${_path_env}")
  endif()
endif()
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the ptscotch header in the given paths
# -------------------------------------------------

set(PAMPA_hdrs_to_find "pampa.h")

# call cmake macro to find the header path
if(PAMPA_INCDIR)
  foreach(pampa_hdr ${PAMPA_hdrs_to_find})
    set(PAMPA_${pampa_hdr}_DIRS "PAMPA_${pampa_hdr}_DIRS-NOTFOUND")
    find_path(PAMPA_${pampa_hdr}_DIRS
      NAMES ${pampa_hdr}
      HINTS ${PAMPA_INCDIR})
    mark_as_advanced(PAMPA_${pampa_hdr}_DIRS)
  endforeach()
else()
  if(PAMPA_DIR)
    foreach(pampa_hdr ${PAMPA_hdrs_to_find})
      set(PAMPA_${pampa_hdr}_DIRS "PAMPA_${pampa_hdr}_DIRS-NOTFOUND")
      find_path(PAMPA_${pampa_hdr}_DIRS
        NAMES ${pampa_hdr}
        HINTS ${PAMPA_DIR}
        PATH_SUFFIXES "include" "include/pampa")
      mark_as_advanced(PAMPA_${pampa_hdr}_DIRS)
    endforeach()
  else()
    foreach(pampa_hdr ${PAMPA_hdrs_to_find})
      set(PAMPA_${pampa_hdr}_DIRS "PAMPA_${pampa_hdr}_DIRS-NOTFOUND")
      find_path(PAMPA_${pampa_hdr}_DIRS
        NAMES ${pampa_hdr}
        HINTS ${_inc_env}
        PATH_SUFFIXES "pampa")
      mark_as_advanced(PAMPA_${pampa_hdr}_DIRS)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
foreach(pampa_hdr ${PAMPA_hdrs_to_find})
  if (PAMPA_${pampa_hdr}_DIRS)
    list(APPEND PAMPA_INCLUDE_DIRS "${PAMPA_${pampa_hdr}_DIRS}")
  else ()
    set(PAMPA_INCLUDE_DIRS "PAMPA_INCLUDE_DIRS-NOTFOUND")
    if (NOT PAMPA_FIND_QUIETLY)
      message(STATUS "Looking for pampa -- ${pampa_hdr} not found")
    endif()
  endif()
endforeach()

# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_PAMPA_LIBDIR "$ENV{PAMPA_LIBDIR}")
if(ENV_PAMPA_LIBDIR)
  list(APPEND _lib_env "${ENV_PAMPA_LIBDIR}")
elseif(ENV_PAMPA_DIR)
  list(APPEND _lib_env "${ENV_PAMPA_DIR}")
  list(APPEND _lib_env "${ENV_PAMPA_DIR}/lib")
else()
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  else()
    if(APPLE)
      string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    else()
      string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    endif()
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  endif()
endif()
list(REMOVE_DUPLICATES _lib_env)

# Try to find the ptscotch lib in the given paths
# ----------------------------------------------

set(PAMPA_libs_to_find "pampa;pampaerr")

# call cmake macro to find the lib path
if(PAMPA_LIBDIR)
  foreach(pampa_lib ${PAMPA_libs_to_find})
    set(PAMPA_${pampa_lib}_LIBRARY "PAMPA_${pampa_lib}_LIBRARY-NOTFOUND")
    find_library(PAMPA_${pampa_lib}_LIBRARY
      NAMES ${pampa_lib}
      HINTS ${PAMPA_LIBDIR})
  endforeach()
else()
  if(PAMPA_DIR)
    foreach(pampa_lib ${PAMPA_libs_to_find})
      set(PAMPA_${pampa_lib}_LIBRARY "PAMPA_${pampa_lib}_LIBRARY-NOTFOUND")
      find_library(PAMPA_${pampa_lib}_LIBRARY
        NAMES ${pampa_lib}
        HINTS ${PAMPA_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    endforeach()
  else()
    foreach(pampa_lib ${PAMPA_libs_to_find})
      set(PAMPA_${pampa_lib}_LIBRARY "PAMPA_${pampa_lib}_LIBRARY-NOTFOUND")
      find_library(PAMPA_${pampa_lib}_LIBRARY
        NAMES ${pampa_lib}
        HINTS ${_lib_env})
    endforeach()
  endif()
endif()

set(PAMPA_LIBRARIES "")
set(PAMPA_LIBRARY_DIRS "")
# If found, add path to cmake variable
# ------------------------------------
foreach(pampa_lib ${PAMPA_libs_to_find})

  if (PAMPA_${pampa_lib}_LIBRARY)
    get_filename_component(${pampa_lib}_lib_path "${PAMPA_${pampa_lib}_LIBRARY}" PATH)
    # set cmake variables
    list(APPEND PAMPA_LIBRARIES "${PAMPA_${pampa_lib}_LIBRARY}")
    list(APPEND PAMPA_LIBRARY_DIRS "${${pampa_lib}_lib_path}")
  else ()
    list(APPEND PAMPA_LIBRARIES "${PAMPA_${pampa_lib}_LIBRARY}")
    if (NOT PAMPA_FIND_QUIETLY)
      message(STATUS "Looking for ptscotch -- lib ${pampa_lib} not found")
    endif()
  endif ()

  mark_as_advanced(PAMPA_${pampa_lib}_LIBRARY)

endforeach()
list(REMOVE_DUPLICATES PAMPA_LIBRARY_DIRS)

# check a function to validate the find
if(PAMPA_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PAMPA
  if (PAMPA_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS  "${PAMPA_INCLUDE_DIRS}")
  endif()
  if (PAMPA_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PAMPA_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PAMPA_LIBRARIES}")
  add_definitions(-Drestrict=__restrict) # pampa uses the restrict keyword
  list(APPEND REQUIRED_FLAGS "-Drestrict=__restrict")
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
    endif()
    if (MPI_C_LINK_FLAGS)
      if (${MPI_C_LINK_FLAGS} MATCHES "  -")
        string(REGEX REPLACE " -" "-" MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS})
      endif()
      list(APPEND REQUIRED_LDFLAGS "${MPI_C_LINK_FLAGS}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # PTSCOTCH
  if (PTSCOTCH_INCLUDE_DIRS_DEP)
    list(APPEND REQUIRED_INCDIRS "${PTSCOTCH_INCLUDE_DIRS_DEP}")
  endif()
  if (PTSCOTCH_CFLAGS_OTHER_DEP)
    list(APPEND REQUIRED_FLAGS "${PTSCOTCH_CFLAGS_OTHER_DEP}")
  endif()
  if (PTSCOTCH_LDFLAGS_OTHER_DEP)
    list(APPEND REQUIRED_LDFLAGS "${PTSCOTCH_LDFLAGS_OTHER_DEP}")
  endif()
  foreach(libdir ${PTSCOTCH_LIBRARY_DIRS_DEP})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  list(APPEND REQUIRED_LIBS "${PTSCOTCH_LIBRARIES_DEP}")
  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  if (REQUIRED_FLAGS)
    set(REQUIRED_FLAGS_COPY "${REQUIRED_FLAGS}")
    set(REQUIRED_FLAGS)
    set(REQUIRED_DEFINITIONS)
    foreach(_flag ${REQUIRED_FLAGS_COPY})
      if (_flag MATCHES "^-D")
       list(APPEND REQUIRED_DEFINITIONS "${_flag}")
      endif()
      string(REGEX REPLACE "^-D.*" "" _flag "${_flag}")
      list(APPEND REQUIRED_FLAGS "${_flag}")
    endforeach()
  endif()
  finds_remove_duplicates()
  set(CMAKE_REQUIRED_DEFINITIONS "${REQUIRED_DEFINITIONS}")
  set(CMAKE_REQUIRED_FLAGS "${REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(PAMPA_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(PAMPA_dmeshInit PAMPA_WORKS)
  mark_as_advanced(PAMPA_WORKS)

  if(PAMPA_WORKS)
    # save link with dependencies
    set(PAMPA_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(PAMPA_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(PAMPA_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(PAMPA_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(PAMPA_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT PAMPA_FIND_QUIETLY)
      message(STATUS "Looking for PAMPA : test of PAMPA_dmeshInit with PAMPA library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(PAMPA_LIBRARIES)

if (PAMPA_LIBRARIES)
  list(GET PAMPA_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT PAMPA_LIBRARY_DIRS)
    set(PAMPA_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PAMPA_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PAMPA library" FORCE)
  else()
    set(PAMPA_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PAMPA library" FORCE)
  endif()
endif()
mark_as_advanced(PAMPA_DIR)
mark_as_advanced(PAMPA_DIR_FOUND)

# Check the size of PAMPA_Num
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${PAMPA_INCLUDE_DIRS})

include(CheckCSourceRuns)
#stdio.h and stdint.h should be included by scotch.h directly
set(PAMPA_C_TEST_PAMPA_Num_4 "
#include <stdio.h>
#include <stdint.h>
#include <pampa.h>
int main(int argc, char **argv) {
  if (sizeof(PAMPA_Num) == 4)
    return 0;
  else
    return 1;
}
")

set(PAMPA_C_TEST_PAMPA_Num_8 "
#include <stdio.h>
#include <stdint.h>
#include <pampa.h>
int main(int argc, char **argv) {
  if (sizeof(PAMPA_Num) == 8)
    return 0;
  else
    return 1;
}
")
check_c_source_runs("${PAMPA_C_TEST_PAMPA_Num_4}" PAMPA_Num_4)
if(NOT PAMPA_Num_4)
  check_c_source_runs("${PAMPA_C_TEST_PAMPA_Num_8}" PAMPA_Num_8)
  if(NOT PAMPA_Num_8)
    set(PAMPA_INTSIZE -1)
  else()
    set(PAMPA_INTSIZE 8)
  endif()
else()
  set(PAMPA_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

# check that PAMPA has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAMPA DEFAULT_MSG
  PAMPA_LIBRARIES
  PAMPA_WORKS)
