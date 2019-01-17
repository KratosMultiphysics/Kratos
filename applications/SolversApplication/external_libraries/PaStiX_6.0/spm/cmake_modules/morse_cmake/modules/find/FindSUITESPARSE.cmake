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
# - Find SUITESPARSE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SUITESPARSE
#               [REQUIRED] # Fail with error if suitesparse is not found
#              )
#
#  SUITESPARSE depends on the following libraries:
#   - BLAS
#   - LAPACK
#   - METIS
##
# This module finds headers and suitesparse library.
# Results are reported in variables:
#  SUITESPARSE_FOUND             - True if headers and requested libraries were found
#  SUITESPARSE_CFLAGS_OTHER      - suitesparse compiler flags without headers paths
#  SUITESPARSE_LDFLAGS_OTHER     - suitesparse linker flags without libraries
#  SUITESPARSE_INCLUDE_DIRS      - suitesparse include directories
#  SUITESPARSE_LIBRARY_DIRS      - suitesparse link directories
#  SUITESPARSE_LIBRARIES         - suitesparse libraries to be linked (absolute path)
#  SUITESPARSE_CFLAGS_OTHER_DEP  - suitesparse + dependencies compiler flags without headers paths
#  SUITESPARSE_LDFLAGS_OTHER_DEP - suitesparse + dependencies linker flags without libraries
#  SUITESPARSE_INCLUDE_DIRS_DEP  - suitesparse + dependencies include directories
#  SUITESPARSE_LIBRARY_DIRS_DEP  - suitesparse + dependencies link directories
#  SUITESPARSE_LIBRARIES_DEP     - suitesparse + dependencies libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DSUITESPARSE_DIR=path/to/suitesparse):
#  SUITESPARSE_DIR              - Where to find the base directory of suitesparse
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SUITESPARSE_DIR

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

if (NOT SUITESPARSE_FOUND)
  set(SUITESPARSE_DIR "" CACHE PATH "Installation directory of SUITESPARSE library")
  if (NOT SUITESPARSE_FIND_QUIETLY)
    message(STATUS "A cache variable, namely SUITESPARSE_DIR, has been set to specify the install directory of SUITESPARSE")
  endif()
endif()

if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE")
endif()

if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - PkgConfig not used")
endif()

# Required dependencies
# ---------------------

if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - Try to detect metis")
endif()
if (SUITESPARSE_FIND_REQUIRED)
  find_package(METIS REQUIRED)
else()
  find_package(METIS)
endif()

# SUITESPARSE depends on LAPACK
#------------------------------
if (NOT SUITESPARSE_FIND_QUIETLY)
  message(STATUS "Looking for SUITESPARSE - Try to detect LAPACK")
endif()
if (SUITESPARSE_FIND_REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()

# Looking for SUITESPARSE
# -----------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_SUITESPARSE_DIR "$ENV{SUITESPARSE_DIR}")
set(ENV_SUITESPARSE_INCDIR "$ENV{SUITESPARSE_INCDIR}")
if(ENV_SUITESPARSE_INCDIR)
  list(APPEND _inc_env "${ENV_SUITESPARSE_INCDIR}")
elseif(ENV_SUITESPARSE_DIR)
  list(APPEND _inc_env "${ENV_SUITESPARSE_DIR}")
  list(APPEND _inc_env "${ENV_SUITESPARSE_DIR}/include")
  list(APPEND _inc_env "${ENV_SUITESPARSE_DIR}/include/suitesparse")
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

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_SUITESPARSE_LIBDIR "$ENV{SUITESPARSE_LIBDIR}")
if(ENV_SUITESPARSE_LIBDIR)
  list(APPEND _lib_env "${ENV_SUITESPARSE_LIBDIR}")
elseif(ENV_SUITESPARSE_DIR)
  list(APPEND _lib_env "${ENV_SUITESPARSE_DIR}")
  list(APPEND _lib_env "${ENV_SUITESPARSE_DIR}/lib")
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

# Looking for include
# -------------------

# Try to find the suitesparse header in the given path
# ----------------------------------------------

# create list of headers to find
list(APPEND SUITESPARSE_hdrs_to_find
  "amd.h"
  "btf.h"
  "ccolamd.h"
  "colamd.h"
  "cs.h"
  "klu.h"
  "ldl.h"
  #"RBio.h"
  "spqr.hpp"
  "SuiteSparse_config.h"
  "umfpack.h")

# call cmake macro to find the header path
if(SUITESPARSE_INCDIR)
  foreach(suitesparse_hdr ${SUITESPARSE_hdrs_to_find})
    set(SUITESPARSE_${suitesparse_hdr}_DIRS "SUITESPARSE_${suitesparse_hdr}_INCLUDE_DIRS-NOTFOUND")
    find_path(SUITESPARSE_${suitesparse_hdr}_DIRS
      NAMES ${suitesparse_hdr}
      HINTS ${SUITESPARSE_INCDIR})
    mark_as_advanced(SUITESPARSE_${suitesparse_hdr}_DIRS)
  endforeach()
else()
  if(SUITESPARSE_DIR)
    set(SUITESPARSE_${suitesparse_hdr}_DIRS "SUITESPARSE_${suitesparse_hdr}_INCLUDE_DIRS-NOTFOUND")
    foreach(suitesparse_hdr ${SUITESPARSE_hdrs_to_find})
      find_path(SUITESPARSE_${suitesparse_hdr}_DIRS
        NAMES ${suitesparse_hdr}
        HINTS ${SUITESPARSE_DIR}
        PATH_SUFFIXES "include")
      mark_as_advanced(SUITESPARSE_${suitesparse_hdr}_DIRS)
    endforeach()
  else()
    foreach(suitesparse_hdr ${SUITESPARSE_hdrs_to_find})
      set(SUITESPARSE_${suitesparse_hdr}_DIRS "SUITESPARSE_${suitesparse_hdr}_INCLUDE_DIRS-NOTFOUND")
      find_path(SUITESPARSE_${suitesparse_hdr}_DIRS
        NAMES ${suitesparse_hdr}
        HINTS ${_inc_env} ${_inc_env}/suitesparse)
      mark_as_advanced(SUITESPARSE_${suitesparse_hdr}_DIRS)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
# detect which precisions are available
foreach(suitesparse_hdr ${SUITESPARSE_hdrs_to_find})
  if (SUITESPARSE_${suitesparse_hdr}_DIRS)
    list(APPEND SUITESPARSE_INCLUDE_DIRS "${SUITESPARSE_${suitesparse_hdr}_DIRS}")
  else ()
    # SuiteSparse_config.h is the minimum to find, consider others as optional?
    if (NOT SUITESPARSE_${SuiteSparse_config.h}_DIRS )
      set(SUITESPARSE_INCLUDE_DIRS "SUITESPARSE_INCLUDE_DIRS-NOTFOUND")
    endif()
    if (NOT SUITESPARSE_FIND_QUIETLY)
      message(STATUS "Looking for suitesparse -- ${suitesparse_hdr} not found")
    endif()
  endif()
endforeach()


# Looking for lib
# ---------------

# create list of libs to find
set(SUITESPARSE_libs_to_find
  "cholmod"
  "cxsparse"
  "klu"
  "ldl"
  "spqr"
  "umfpack"
  "amd"
  "btf"
  "camd"
  "ccolamd"
  "colamd"
  #"rbio"
  "suitesparseconfig"
  )

# call cmake macro to find the lib path
if(SUITESPARSE_LIBDIR)
  foreach(suitesparse_lib ${SUITESPARSE_libs_to_find})
    set(SUITESPARSE_${suitesparse_lib}_LIBRARY "SUITESPARSE_${suitesparse_lib}_LIBRARY-NOTFOUND")
    find_library(SUITESPARSE_${suitesparse_lib}_LIBRARY
      NAMES ${suitesparse_lib}
      HINTS ${SUITESPARSE_LIBDIR})
    mark_as_advanced(SUITESPARSE_${suitesparse_lib}_LIBRARY)
  endforeach()
else()
  if(SUITESPARSE_DIR)
    foreach(suitesparse_lib ${SUITESPARSE_libs_to_find})
      set(SUITESPARSE_${suitesparse_lib}_LIBRARY "SUITESPARSE_${suitesparse_lib}_LIBRARY-NOTFOUND")
      find_library(SUITESPARSE_${suitesparse_lib}_LIBRARY
        NAMES ${suitesparse_lib}
        HINTS ${SUITESPARSE_DIR}
        PATH_SUFFIXES lib lib32 lib64)
      mark_as_advanced(SUITESPARSE_${suitesparse_lib}_LIBRARY)
    endforeach()
  else()
    foreach(suitesparse_lib ${SUITESPARSE_libs_to_find})
      set(SUITESPARSE_${suitesparse_lib}_LIBRARY "SUITESPARSE_${suitesparse_lib}_LIBRARY-NOTFOUND")
      find_library(SUITESPARSE_${suitesparse_lib}_LIBRARY
        NAMES ${suitesparse_lib}
        HINTS ${_lib_env})
      mark_as_advanced(SUITESPARSE_${suitesparse_lib}_LIBRARY)
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
set(SUITESPARSE_LIBRARIES "")
set(SUITESPARSE_LIBRARY_DIRS "")
foreach(suitesparse_lib ${SUITESPARSE_libs_to_find})
  if (SUITESPARSE_${suitesparse_lib}_LIBRARY)
    list(APPEND SUITESPARSE_LIBRARIES "${SUITESPARSE_${suitesparse_lib}_LIBRARY}")
    get_filename_component(${suitesparse_lib}_lib_path ${SUITESPARSE_${suitesparse_lib}_LIBRARY} PATH)
    list(APPEND SUITESPARSE_LIBRARY_DIRS "${${suitesparse_lib}_lib_path}")
  else ()
    # libsuitesparseconfig is the minimum to find, consider others as optional?
    if (NOT SUITESPARSE_suitesparseconfig_LIBRARY)
      list(APPEND SUITESPARSE_LIBRARIES "${SUITESPARSE_LIBRARIES-NOTFOUND}")
    endif()
    if(NOT SUITESPARSE_FIND_QUIETLY)
      message(STATUS "Looking for suitesparse -- lib${suitesparse_lib}.a/so not found")
    endif()
  endif()
endforeach()
list(REMOVE_DUPLICATES SUITESPARSE_LIBRARY_DIRS)
list(REMOVE_DUPLICATES SUITESPARSE_INCLUDE_DIRS)

# check a function to validate the find
if(SUITESPARSE_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # SUITESPARSE
  if (SUITESPARSE_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${SUITESPARSE_INCLUDE_DIRS}")
  endif()
  foreach(libdir ${SUITESPARSE_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${SUITESPARSE_LIBRARIES}")
  # METIS
  if (METIS_FOUND)
    if (METIS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${METIS_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${METIS_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${METIS_LIBRARIES}")
  endif()
  # LAPACK
  if (LAPACK_FOUND)
    if (LAPACK_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS_DEP}")
    endif()
    if (LAPACK_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${LAPACK_CFLAGS_OTHER_DEP}")
    endif()
    if (LAPACK_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${LAPACK_LDFLAGS_OTHER_DEP}")
    endif()
    if (LAPACK_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS_DEP}")
    endif()
    list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES_DEP}")
  endif()
  # others
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
  endif()

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
  include(CheckFortranFunctionExists)
  unset(SUITESPARSE_WORKS CACHE)
  check_function_exists(SuiteSparse_start SUITESPARSE_WORKS)
  mark_as_advanced(SUITESPARSE_WORKS)

  if(SUITESPARSE_WORKS)
    # save link with dependencies
    set(SUITESPARSE_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(SUITESPARSE_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(SUITESPARSE_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(SUITESPARSE_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(SUITESPARSE_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT SUITESPARSE_FIND_QUIETLY)
      message(STATUS "Looking for SUITESPARSE : test of symbol SuiteSparse_start fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe SUITESPARSE is linked with specific libraries. ")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(SUITESPARSE_LIBRARIES)

if (SUITESPARSE_LIBRARIES)
  list(GET SUITESPARSE_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT SUITESPARSE_LIBRARY_DIRS)
    set(SUITESPARSE_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(SUITESPARSE_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of SUITESPARSE library" FORCE)
  else()
    set(SUITESPARSE_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of SUITESPARSE library" FORCE)
  endif()
endif()
mark_as_advanced(SUITESPARSE_DIR)
mark_as_advanced(SUITESPARSE_DIR_FOUND)

# check that SUITESPARSE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUITESPARSE DEFAULT_MSG
  SUITESPARSE_LIBRARIES
  SUITESPARSE_WORKS)
