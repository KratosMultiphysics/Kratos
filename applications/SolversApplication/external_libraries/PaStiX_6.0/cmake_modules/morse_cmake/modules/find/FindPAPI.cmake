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
# - Find PAPI include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PAPI
#               [REQUIRED]) # Fail with error if papi is not found
#
# This module finds headers and papi library.
# Results are reported in variables:
#  PAPI_FOUND           - True if headers and requested libraries were found
#  PAPI_INCLUDE_DIRS    - papi include directories
#  PAPI_LIBRARY_DIRS    - Link directories for papi libraries
#  PAPI_LIBRARIES       - papi component libraries to be linked
#
#  PAPI_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = PAPI
#  <XPREFIX> = <PREFIX>        for common case
#  <XPREFIX> = <PREFIX>_STATIC for static linking
#  <XPREFIX>_FOUND          ... set to 1 if module(s) exist
#  <XPREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
#  <XPREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
#  <XPREFIX>_LDFLAGS        ... all required linker flags
#  <XPREFIX>_LDFLAGS_OTHER  ... all other linker flags
#  <XPREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
#  <XPREFIX>_CFLAGS         ... all required cflags
#  <XPREFIX>_CFLAGS_OTHER   ... the other compiler flags
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DPAPI_DIR=path/to/papi):
#  PAPI_DIR             - Where to find the base directory of papi
#  PAPI_INCDIR          - Where to find the header files
#  PAPI_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PAPI_DIR, PAPI_INCDIR, PAPI_LIBDIR

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

if (NOT PAPI_FOUND)
  set(PAPI_DIR "" CACHE PATH "Installation directory of PAPI library")
  if (NOT PAPI_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PAPI_DIR, has been set to specify the install directory of PAPI")
  endif()
endif()

set(ENV_PAPI_DIR "$ENV{PAPI_DIR}")
set(ENV_PAPI_INCDIR "$ENV{PAPI_INCDIR}")
set(ENV_PAPI_LIBDIR "$ENV{PAPI_LIBDIR}")
set(PAPI_GIVEN_BY_USER "FALSE")
if ( PAPI_DIR OR ( PAPI_INCDIR AND PAPI_LIBDIR) OR ENV_PAPI_DIR OR (ENV_PAPI_INCDIR AND ENV_PAPI_LIBDIR) )
  set(PAPI_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT PAPI_GIVEN_BY_USER)

  pkg_search_module(PAPI papi)

  if (NOT PAPI_FIND_QUIETLY)
    if (PAPI_FOUND AND PAPI_LIBRARIES)
      message(STATUS "Looking for PAPI - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for PAPI - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing papi.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (PAPI_FOUND AND PAPI_LIBRARIES)
    set(PAPI_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(PAPI)
  else()
    set(PAPI_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT PAPI_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PAPI_FOUND) OR (PAPI_GIVEN_BY_USER) )

  if (NOT PAPI_FIND_QUIETLY)
    message(STATUS "Looking for PAPI - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_PAPI_DIR "$ENV{PAPI_DIR}")
  set(ENV_PAPI_INCDIR "$ENV{PAPI_INCDIR}")
  if(ENV_PAPI_INCDIR)
    list(APPEND _inc_env "${ENV_PAPI_INCDIR}")
  elseif(ENV_PAPI_DIR)
    list(APPEND _inc_env "${ENV_PAPI_DIR}")
    list(APPEND _inc_env "${ENV_PAPI_DIR}/include")
    list(APPEND _inc_env "${ENV_PAPI_DIR}/include/papi")
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

  # Try to find the papi header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(PAPI_INCDIR)
    set(PAPI_papi.h_DIRS "PAPI_papi.h_DIRS-NOTFOUND")
    find_path(PAPI_papi.h_DIRS
      NAMES papi.h
      HINTS ${PAPI_INCDIR})
  else()
    if(PAPI_DIR)
      set(PAPI_papi.h_DIRS "PAPI_papi.h_DIRS-NOTFOUND")
      find_path(PAPI_papi.h_DIRS
        NAMES papi.h
        HINTS ${PAPI_DIR}
        PATH_SUFFIXES "include" "include/papi")
    else()
      set(PAPI_papi.h_DIRS "PAPI_papi.h_DIRS-NOTFOUND")
      find_path(PAPI_papi.h_DIRS
        NAMES papi.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "papi")
    endif()
  endif()
  mark_as_advanced(PAPI_papi.h_DIRS)

  # Add path to cmake variable
  # ------------------------------------
  if (PAPI_papi.h_DIRS)
    set(PAPI_INCLUDE_DIRS "${PAPI_papi.h_DIRS}")
  else ()
    set(PAPI_INCLUDE_DIRS "PAPI_INCLUDE_DIRS-NOTFOUND")
    if(NOT PAPI_FIND_QUIETLY)
      message(STATUS "Looking for papi -- papi.h not found")
    endif()
  endif ()

  if (PAPI_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES PAPI_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_PAPI_LIBDIR "$ENV{PAPI_LIBDIR}")
  if(ENV_PAPI_LIBDIR)
    list(APPEND _lib_env "${ENV_PAPI_LIBDIR}")
  elseif(ENV_PAPI_DIR)
    list(APPEND _lib_env "${ENV_PAPI_DIR}")
    list(APPEND _lib_env "${ENV_PAPI_DIR}/lib")
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

  # Try to find the papi lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(PAPI_LIBDIR)
    set(PAPI_papi_LIBRARY "PAPI_papi_LIBRARY-NOTFOUND")
    find_library(PAPI_papi_LIBRARY
      NAMES papi
      HINTS ${PAPI_LIBDIR})
  else()
    if(PAPI_DIR)
      set(PAPI_papi_LIBRARY "PAPI_papi_LIBRARY-NOTFOUND")
      find_library(PAPI_papi_LIBRARY
        NAMES papi
        HINTS ${PAPI_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(PAPI_papi_LIBRARY "PAPI_papi_LIBRARY-NOTFOUND")
      find_library(PAPI_papi_LIBRARY
        NAMES papi
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(PAPI_papi_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (PAPI_papi_LIBRARY)
    get_filename_component(papi_lib_path ${PAPI_papi_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(PAPI_LIBRARIES    "${PAPI_papi_LIBRARY}")
    set(PAPI_LIBRARY_DIRS "${papi_lib_path}")
  else ()
    set(PAPI_LIBRARIES    "PAPI_LIBRARIES-NOTFOUND")
    set(PAPI_LIBRARY_DIRS "PAPI_LIBRARY_DIRS-NOTFOUND")
    if(NOT PAPI_FIND_QUIETLY)
      message(STATUS "Looking for papi -- lib papi not found")
    endif()
  endif ()

  if (PAPI_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES PAPI_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PAPI_FOUND) OR (PAPI_GIVEN_BY_USER) )

# check a function to validate the find
if(PAPI_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PAPI
  if (PAPI_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${PAPI_INCLUDE_DIRS}")
  endif()
  if (PAPI_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${PAPI_CFLAGS_OTHER}")
  endif()
  if (PAPI_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${PAPI_LDFLAGS_OTHER}")
  endif()
  if (PAPI_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${PAPI_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${PAPI_LIBRARIES}")

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
  unset(PAPI_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(PAPI_start PAPI_WORKS)
  mark_as_advanced(PAPI_WORKS)

  if(NOT PAPI_WORKS)
    if(NOT PAPI_FIND_QUIETLY)
      message(STATUS "Looking for papi : test of PAPI_start with papi library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(PAPI_LIBRARIES)

if (PAPI_LIBRARIES)
  if (PAPI_LIBRARY_DIRS)
    list(GET PAPI_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET PAPI_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(PAPI_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PAPI_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PAPI library" FORCE)
  else()
    set(PAPI_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PAPI library" FORCE)
  endif()
endif()
mark_as_advanced(PAPI_DIR)
mark_as_advanced(PAPI_DIR_FOUND)

# check that PAPI has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PAPI DEFAULT_MSG
  PAPI_LIBRARIES
  PAPI_WORKS)
