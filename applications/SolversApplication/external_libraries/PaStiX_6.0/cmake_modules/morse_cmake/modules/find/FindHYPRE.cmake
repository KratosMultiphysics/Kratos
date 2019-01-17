###
# WARNING: only HYPRE lib is searched for now
# it is surely too simple, must be completed
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
# - Find HYPRE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HYPRE
#               [REQUIRED]) # Fail with error if hypre is not found
#
#  HYPRE depends on the following libraries:
#   - MPI
#
# This module finds headers and hypre library.
# Results are reported in variables:
#  HYPRE_FOUND             - True if headers and requested libraries were found
#  HYPRE_INCLUDE_DIRS      - hypre include directories
#  HYPRE_LIBRARY_DIRS      - Link directories for hypre libraries
#  HYPRE_LIBRARIES         - hypre component libraries to be linked
#  HYPRE_INCLUDE_DIRS_DEP  - hypre + dependencies include directories
#  HYPRE_LIBRARY_DIRS_DEP  - hypre + dependencies link directories
#  HYPRE_LIBRARIES_DEP     - hypre + dependencies libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DHYPRE_DIR=path/to/hypre):
#  HYPRE_DIR             - Where to find the base directory of hypre
#  HYPRE_INCDIR          - Where to find the header files
#  HYPRE_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: HYPRE_DIR, HYPRE_INCDIR, HYPRE_LIBDIR

#=============================================================================
# Copyright 2016-2018 Inria
# Copyright 2016-2018 Florent Pruvost
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

if (NOT HYPRE_FOUND)
  set(HYPRE_DIR "" CACHE PATH "Installation directory of HYPRE library")
  if (NOT HYPRE_FIND_QUIETLY)
    message(STATUS "A cache variable, namely HYPRE_DIR, has been set to specify the install directory of HYPRE")
  endif()
endif()

set(ENV_HYPRE_DIR "$ENV{HYPRE_DIR}")
set(ENV_HYPRE_INCDIR "$ENV{HYPRE_INCDIR}")
set(ENV_HYPRE_LIBDIR "$ENV{HYPRE_LIBDIR}")
set(HYPRE_GIVEN_BY_USER "FALSE")
if ( HYPRE_DIR OR ( HYPRE_INCDIR AND HYPRE_LIBDIR) OR ENV_HYPRE_DIR OR (ENV_HYPRE_INCDIR AND ENV_HYPRE_LIBDIR) )
  set(HYPRE_GIVEN_BY_USER "TRUE")
endif()

# HYPRE depends on MPI, try to find it
if (HYPRE_FIND_REQUIRED)
  find_package(MPI REQUIRED)
else()
  find_package(MPI)
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
if(ENV_HYPRE_INCDIR)
  list(APPEND _inc_env "${ENV_HYPRE_INCDIR}")
elseif(ENV_HYPRE_DIR)
  list(APPEND _inc_env "${ENV_HYPRE_DIR}")
  list(APPEND _inc_env "${ENV_HYPRE_DIR}/include")
  list(APPEND _inc_env "${ENV_HYPRE_DIR}/include/hypre")
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

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_inc_env}")

# Try to find the hypre header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(HYPRE_INCDIR)
  set(HYPRE_HYPRE.h_DIRS "HYPRE_HYPRE.h_DIRS-NOTFOUND")
  find_path(HYPRE_HYPRE.h_DIRS
    NAMES HYPRE.h
    HINTS ${HYPRE_INCDIR})
else()
  if(HYPRE_DIR)
    set(HYPRE_HYPRE.h_DIRS "HYPRE_HYPRE.h_DIRS-NOTFOUND")
    find_path(HYPRE_HYPRE.h_DIRS
      NAMES HYPRE.h
      HINTS ${HYPRE_DIR}
      PATH_SUFFIXES "include" "include/hypre")
  else()
    set(HYPRE_HYPRE.h_DIRS "HYPRE_HYPRE.h_DIRS-NOTFOUND")
    find_path(HYPRE_HYPRE.h_DIRS
      NAMES HYPRE.h
      HINTS ${PATH_TO_LOOK_FOR}
      PATH_SUFFIXES "hypre")
  endif()
endif()
mark_as_advanced(HYPRE_HYPRE.h_DIRS)

# Add path to cmake variable
# ------------------------------------
if (HYPRE_HYPRE.h_DIRS)
  set(HYPRE_INCLUDE_DIRS "${HYPRE_HYPRE.h_DIRS}")
else ()
  set(HYPRE_INCLUDE_DIRS "HYPRE_INCLUDE_DIRS-NOTFOUND")
  if(NOT HYPRE_FIND_QUIETLY)
    message(STATUS "Looking for hypre -- HYPRE.h not found")
  endif()
endif ()

if (HYPRE_INCLUDE_DIRS)
  list(REMOVE_DUPLICATES HYPRE_INCLUDE_DIRS)
endif ()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
if(ENV_HYPRE_LIBDIR)
  list(APPEND _lib_env "${ENV_HYPRE_LIBDIR}")
elseif(ENV_HYPRE_DIR)
  list(APPEND _lib_env "${ENV_HYPRE_DIR}")
  list(APPEND _lib_env "${ENV_HYPRE_DIR}/lib")
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

# set paths where to look for
set(PATH_TO_LOOK_FOR "${_lib_env}")

# Try to find the hypre lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(HYPRE_LIBDIR)
  set(HYPRE_HYPRE_LIBRARY "HYPRE_HYPRE_LIBRARY-NOTFOUND")
  find_library(HYPRE_HYPRE_LIBRARY
    NAMES HYPRE
    HINTS ${HYPRE_LIBDIR})
else()
  if(HYPRE_DIR)
    set(HYPRE_HYPRE_LIBRARY "HYPRE_HYPRE_LIBRARY-NOTFOUND")
    find_library(HYPRE_HYPRE_LIBRARY
      NAMES HYPRE
      HINTS ${HYPRE_DIR}
      PATH_SUFFIXES lib lib32 lib64)
  else()
    set(HYPRE_HYPRE_LIBRARY "HYPRE_HYPRE_LIBRARY-NOTFOUND")
    find_library(HYPRE_HYPRE_LIBRARY
      NAMES HYPRE
      HINTS ${PATH_TO_LOOK_FOR})
  endif()
endif()
mark_as_advanced(HYPRE_HYPRE_LIBRARY)

# If found, add path to cmake variable
# ------------------------------------
if (HYPRE_HYPRE_LIBRARY)
  get_filename_component(hypre_lib_path ${HYPRE_HYPRE_LIBRARY} PATH)
  # set cmake variables (respects naming convention)
  set(HYPRE_LIBRARIES    "${HYPRE_HYPRE_LIBRARY}")
  set(HYPRE_LIBRARY_DIRS "${hypre_lib_path}")
else ()
  set(HYPRE_LIBRARIES    "HYPRE_LIBRARIES-NOTFOUND")
  set(HYPRE_LIBRARY_DIRS "HYPRE_LIBRARY_DIRS-NOTFOUND")
  if(NOT HYPRE_FIND_QUIETLY)
    message(STATUS "Looking for hypre -- lib HYPRE not found")
  endif()
endif ()

if (HYPRE_LIBRARY_DIRS)
  list(REMOVE_DUPLICATES HYPRE_LIBRARY_DIRS)
endif ()

# check a function to validate the find
if(HYPRE_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # HYPRE
  if (HYPRE_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${HYPRE_INCLUDE_DIRS}")
  endif()
  if (HYPRE_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${HYPRE_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${HYPRE_LIBRARIES}")
  # MPI
  if (MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    if (MPI_C_LINK_FLAGS)
      if (${MPI_C_LINK_FLAGS} MATCHES "  -")
        string(REGEX REPLACE " -" "-" MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS})
      endif()
      list(APPEND REQUIRED_LDFLAGS "${MPI_C_LINK_FLAGS}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # libm
  set(M_LIBRARY "M_LIBRARY-NOTFOUND")
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "${M_LIBRARY}")
  endif()
  # libstdc++
  set(stdcpp_LIBRARY "stdcpp_LIBRARY-NOTFOUND")
  find_library(stdcpp_LIBRARY NAMES stdc++)
  mark_as_advanced(stdcpp_LIBRARY)
  if(stdcpp_LIBRARY)
    list(APPEND REQUIRED_LIBS "${stdcpp_LIBRARY}")
  endif()

  # set required libraries for link
  finds_remove_duplicates()
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(HYPRE_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(HYPRE_StructGridCreate HYPRE_WORKS)
  mark_as_advanced(HYPRE_WORKS)

  if(HYPRE_WORKS)
    # save link with dependencies
    set(HYPRE_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(HYPRE_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(HYPRE_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
  else()
    if(NOT HYPRE_FIND_QUIETLY)
      message(STATUS "Looking for HYPRE : test of HYPRE_StructGridCreate with HYPRE library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(HYPRE_LIBRARIES)

if (HYPRE_LIBRARIES)
  if (HYPRE_LIBRARY_DIRS)
    list(GET HYPRE_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET HYPRE_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(HYPRE_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of HYPRE library" FORCE)
  else()
    set(HYPRE_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of HYPRE library" FORCE)
  endif()
endif()
mark_as_advanced(HYPRE_DIR)
mark_as_advanced(HYPRE_DIR_FOUND)

# check that HYPRE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HYPRE DEFAULT_MSG
  HYPRE_LIBRARIES
  HYPRE_WORKS)
