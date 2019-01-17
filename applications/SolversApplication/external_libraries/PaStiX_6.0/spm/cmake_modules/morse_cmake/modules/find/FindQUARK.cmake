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
# - Find QUARK include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(QUARK
#               [REQUIRED]             # Fail with error if quark is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  QUARK depends on the following libraries:
#   - Threads
#
#  COMPONENTS are optional libraries QUARK could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - HWLOC: to activate detection of QUARK linked with hwloc
#
# This module finds headers and quark library.
# Results are reported in variables:
#  QUARK_FOUND             - True if headers and requested libraries were found
#  QUARK_CFLAGS_OTHER      - quark compiler flags without headers paths
#  QUARK_LDFLAGS_OTHER     - quark linker flags without libraries
#  QUARK_INCLUDE_DIRS      - quark include directories
#  QUARK_LIBRARY_DIRS      - quark link directories
#  QUARK_LIBRARIES         - quark libraries to be linked (absolute path)
#  QUARK_CFLAGS_OTHER_DEP  - quark + dependencies compiler flags without headers paths
#  QUARK_LDFLAGS_OTHER_DEP - quark + dependencies linker flags without libraries
#  QUARK_INCLUDE_DIRS_DEP  - quark + dependencies include directories
#  QUARK_LIBRARY_DIRS_DEP  - quark + dependencies link directories
#  QUARK_LIBRARIES_DEP     - quark + dependencies libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DQUARK=path/to/quark):
#  QUARK_DIR              - Where to find the base directory of quark
#  QUARK_INCDIR           - Where to find the header files
#  QUARK_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: QUARK_DIR, QUARK_INCDIR, QUARK_LIBDIR

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

if (NOT QUARK_FOUND)
  set(QUARK_DIR "" CACHE PATH "Installation directory of QUARK library")
  if (NOT QUARK_FIND_QUIETLY)
    message(STATUS "A cache variable, namely QUARK_DIR, has been set to specify the install directory of QUARK")
  endif()
endif()

# QUARK may depend on HWLOC
# try to find it specified as COMPONENTS during the call
set(QUARK_LOOK_FOR_HWLOC FALSE)

if( QUARK_FIND_COMPONENTS )
  foreach( component ${QUARK_FIND_COMPONENTS} )
    if(${component} STREQUAL "HWLOC")
      set(QUARK_LOOK_FOR_HWLOC TRUE)
    endif()
  endforeach()
endif()

# QUARK may depend on Threads, try to find it
if (QUARK_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()

# QUARK may depend on HWLOC, try to find it
if (QUARK_LOOK_FOR_HWLOC)
  if (QUARK_FIND_REQUIRED AND QUARK_FIND_REQUIRED_HWLOC)
    find_package(HWLOC REQUIRED)
  else()
    find_package(HWLOC)
  endif()
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_QUARK_DIR "$ENV{QUARK_DIR}")
set(ENV_QUARK_INCDIR "$ENV{QUARK_INCDIR}")
if(ENV_QUARK_INCDIR)
  list(APPEND _inc_env "${ENV_QUARK_INCDIR}")
elseif(ENV_QUARK_DIR)
  list(APPEND _inc_env "${ENV_QUARK_DIR}")
  list(APPEND _inc_env "${ENV_QUARK_DIR}/include")
  list(APPEND _inc_env "${ENV_QUARK_DIR}/include/quark")
  list(APPEND _inc_env "${ENV_QUARK_DIR}/include/plasma")
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


# Try to find the quark header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(QUARK_INCDIR)
  set(QUARK_quark.h_DIRS "QUARK_quark.h_DIRS-NOTFOUND")
  find_path(QUARK_quark.h_DIRS
    NAMES quark.h
    HINTS ${QUARK_INCDIR})
else()
  if(QUARK_DIR)
    set(QUARK_quark.h_DIRS "QUARK_quark.h_DIRS-NOTFOUND")
    find_path(QUARK_quark.h_DIRS
      NAMES quark.h
      HINTS ${QUARK_DIR}
      PATH_SUFFIXES "include" "include/quark" "include/plasma")
  else()
    set(QUARK_quark.h_DIRS "QUARK_quark.h_DIRS-NOTFOUND")
    find_path(QUARK_quark.h_DIRS
      NAMES quark.h
      HINTS ${_inc_env})
  endif()
endif()
mark_as_advanced(QUARK_quark.h_DIRS)

# If found, add path to cmake variable
# ------------------------------------
if (QUARK_quark.h_DIRS)
  set(QUARK_INCLUDE_DIRS "${QUARK_quark.h_DIRS}")
else ()
  set(QUARK_INCLUDE_DIRS "QUARK_INCLUDE_DIRS-NOTFOUND")
  if(NOT QUARK_FIND_QUIETLY)
    message(STATUS "Looking for quark -- quark.h not found")
  endif()
endif()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_QUARK_LIBDIR "$ENV{QUARK_LIBDIR}")
if(ENV_QUARK_LIBDIR)
  list(APPEND _lib_env "${ENV_QUARK_LIBDIR}")
elseif(ENV_QUARK_DIR)
  list(APPEND _lib_env "${ENV_QUARK_DIR}")
  list(APPEND _lib_env "${ENV_QUARK_DIR}/lib")
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

# Try to find the quark lib in the given paths
# ----------------------------------------------

# call cmake macro to find the lib path
if(QUARK_LIBDIR)
  set(QUARK_quark_LIBRARY "QUARK_quark_LIBRARY-NOTFOUND")
  find_library(QUARK_quark_LIBRARY
    NAMES quark
    HINTS ${QUARK_LIBDIR})
else()
  if(QUARK_DIR)
    set(QUARK_quark_LIBRARY "QUARK_quark_LIBRARY-NOTFOUND")
    find_library(QUARK_quark_LIBRARY
      NAMES quark
      HINTS ${QUARK_DIR}
      PATH_SUFFIXES lib lib32 lib64)
  else()
    set(QUARK_quark_LIBRARY "QUARK_quark_LIBRARY-NOTFOUND")
    find_library(QUARK_quark_LIBRARY
      NAMES quark
      HINTS ${_lib_env})
  endif()
endif()
mark_as_advanced(QUARK_quark_LIBRARY)

# If found, add path to cmake variable
# ------------------------------------
if (QUARK_quark_LIBRARY)
  get_filename_component(quark_lib_path "${QUARK_quark_LIBRARY}" PATH)
  # set cmake variables
  set(QUARK_LIBRARIES    "${QUARK_quark_LIBRARY}")
  set(QUARK_LIBRARY_DIRS "${quark_lib_path}")
else ()
  set(QUARK_LIBRARIES    "QUARK_LIBRARIES-NOTFOUND")
  set(QUARK_LIBRARY_DIRS "QUARK_LIBRARY_DIRS-NOTFOUND")
  if(NOT QUARK_FIND_QUIETLY)
    message(STATUS "Looking for quark -- lib quark not found")
  endif()
endif ()

# check a function to validate the find
if(QUARK_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # QUARK
  if (QUARK_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${QUARK_INCLUDE_DIRS}")
  endif()
  if (QUARK_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${QUARK_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${QUARK_LIBRARIES}")
  # HWLOC
  if (HWLOC_FOUND AND QUARK_LOOK_FOR_HWLOC)
    if (HWLOC_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
    endif()
    if (HWLOC_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${HWLOC_CFLAGS_OTHER}")
    endif()
    if (HWLOC_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${HWLOC_LDFLAGS_OTHER}")
    endif()
    if (HWLOC_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${HWLOC_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${HWLOC_LIBRARIES}")
  endif()
  # THREADS
  list(APPEND REQUIRED_LIBS "${CMAKE_THREAD_LIBS_INIT}")

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
  unset(QUARK_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(QUARK_New QUARK_WORKS)
  mark_as_advanced(QUARK_WORKS)

  if(QUARK_WORKS)
    # save link with dependencies
    set(QUARK_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(QUARK_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(QUARK_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(QUARK_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(QUARK_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT QUARK_FIND_QUIETLY)
      message(STATUS "Looking for QUARK : test of QUARK_New with QUARK library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe QUARK is linked with specific libraries like. "
        "Have you tried with COMPONENTS (HWLOC)? See the explanation in FindQUARK.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(QUARK_LIBRARIES)

if (QUARK_LIBRARIES)
  list(GET QUARK_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT QUARK_LIBRARY_DIRS)
    set(QUARK_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(QUARK_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of QUARK library" FORCE)
  else()
    set(QUARK_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of QUARK library" FORCE)
  endif()
endif()
mark_as_advanced(QUARK_DIR)
mark_as_advanced(QUARK_DIR_FOUND)

# check that QUARK has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QUARK DEFAULT_MSG
  QUARK_LIBRARIES
  QUARK_WORKS)
