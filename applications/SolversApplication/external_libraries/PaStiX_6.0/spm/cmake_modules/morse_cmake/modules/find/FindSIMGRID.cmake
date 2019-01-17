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
# - Find SIMGRID include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SIMGRID
#               [REQUIRED]) # Fail with error if simgrid is not found
#
# This module finds headers and simgrid library.
# Results are reported in variables:
#  SIMGRID_FOUND           - True if headers and requested libraries were found
#  SIMGRID_INCLUDE_DIRS    - simgrid include directories
#  SIMGRID_LIBRARY_DIRS    - Link directories for simgrid libraries
#  SIMGRID_LIBRARIES       - simgrid component libraries to be linked
#
#  SIMGRID_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = SIMGRID
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
# options at configure (ex: cmake path/to/project -DSIMGRID_DIR=path/to/simgrid):
#  SIMGRID_DIR             - Where to find the base directory of simgrid
#  SIMGRID_INCDIR          - Where to find the header files
#  SIMGRID_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SIMGRID_DIR, SIMGRID_INCDIR, SIMGRID_LIBDIR

#=============================================================================
# Copyright 2012-2013 Inria
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

if (NOT SIMGRID_FOUND)
  set(SIMGRID_DIR "" CACHE PATH "Installation directory of SIMGRID library")
  if (NOT SIMGRID_FIND_QUIETLY)
    message(STATUS "A cache variable, namely SIMGRID_DIR, has been set to specify the install directory of SIMGRID")
  endif()
endif()

set(ENV_SIMGRID_DIR "$ENV{SIMGRID_DIR}")
set(ENV_SIMGRID_INCDIR "$ENV{SIMGRID_INCDIR}")
set(ENV_SIMGRID_LIBDIR "$ENV{SIMGRID_LIBDIR}")
set(SIMGRID_GIVEN_BY_USER "FALSE")
if ( SIMGRID_DIR OR ( SIMGRID_INCDIR AND SIMGRID_LIBDIR) OR ENV_SIMGRID_DIR OR (ENV_SIMGRID_INCDIR AND ENV_SIMGRID_LIBDIR) )
  set(SIMGRID_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT SIMGRID_GIVEN_BY_USER)

  pkg_search_module(SIMGRID simgrid)

  if (NOT SIMGRID_FIND_QUIETLY)
    if (SIMGRID_FOUND AND SIMGRID_LIBRARIES)
      message(STATUS "Looking for SIMGRID - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for SIMGRID - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing simgrid.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (SIMGRID_FOUND AND SIMGRID_LIBRARIES)
    set(SIMGRID_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(SIMGRID)
  else()
    set(SIMGRID_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT SIMGRID_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT SIMGRID_FOUND) OR (SIMGRID_GIVEN_BY_USER) )

  if (NOT SIMGRID_FIND_QUIETLY)
    message(STATUS "Looking for SIMGRID - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_SIMGRID_DIR "$ENV{SIMGRID_DIR}")
  set(ENV_SIMGRID_INCDIR "$ENV{SIMGRID_INCDIR}")
  if(ENV_SIMGRID_INCDIR)
    list(APPEND _inc_env "${ENV_SIMGRID_INCDIR}")
  elseif(ENV_SIMGRID_DIR)
    list(APPEND _inc_env "${ENV_SIMGRID_DIR}")
    list(APPEND _inc_env "${ENV_SIMGRID_DIR}/include")
    list(APPEND _inc_env "${ENV_SIMGRID_DIR}/include/simgrid")
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

  # Try to find the simgrid header in the given paths
  # -------------------------------------------------
  set(SIMGRID_hdrs_to_find "simgrid_config.h")

  # call cmake macro to find the header path
  if(SIMGRID_INCDIR)
    foreach(simgrid_hdr ${SIMGRID_hdrs_to_find})
      set(SIMGRID_${simgrid_hdr}_DIRS "SIMGRID_${simgrid_hdr}_DIRS-NOTFOUND")
      find_path(SIMGRID_${simgrid_hdr}_DIRS
        NAMES ${simgrid_hdr}
        HINTS ${SIMGRID_INCDIR})
        mark_as_advanced(SIMGRID_${simgrid_hdr}_DIRS)
    endforeach()
  else()
    if(SIMGRID_DIR)
      foreach(simgrid_hdr ${SIMGRID_hdrs_to_find})
        set(SIMGRID_${simgrid_hdr}_DIRS "SIMGRID_${simgrid_hdr}_DIRS-NOTFOUND")
        find_path(SIMGRID_${simgrid_hdr}_DIRS
          NAMES ${simgrid_hdr}
          HINTS ${SIMGRID_DIR}
          PATH_SUFFIXES "include" "include/simgrid")
          mark_as_advanced(SIMGRID_${simgrid_hdr}_DIRS)
      endforeach()
    else()
      foreach(simgrid_hdr ${SIMGRID_hdrs_to_find})
        set(SIMGRID_${simgrid_hdr}_DIRS "SIMGRID_${simgrid_hdr}_DIRS-NOTFOUND")
        find_path(SIMGRID_${simgrid_hdr}_DIRS
          NAMES ${simgrid_hdr}
          HINTS ${_inc_env}
          PATH_SUFFIXES "simgrid")
          mark_as_advanced(SIMGRID_${simgrid_hdr}_DIRS)
      endforeach()
    endif()
  endif()

  # Add path to cmake variable
  # ------------------------------------
  if (SIMGRID_simgrid_config.h_DIRS)
    set(SIMGRID_INCLUDE_DIRS "${SIMGRID_simgrid_config.h_DIRS}")
  else ()
    set(SIMGRID_INCLUDE_DIRS "SIMGRID_INCLUDE_DIRS-NOTFOUND")
    if(NOT SIMGRID_FIND_QUIETLY)
      message(STATUS "Looking for simgrid -- simgrid_config.h not found")
    endif()
  endif ()

  if (SIMGRID_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES SIMGRID_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_SIMGRID_LIBDIR "$ENV{SIMGRID_LIBDIR}")
  if(ENV_SIMGRID_LIBDIR)
    list(APPEND _lib_env "${ENV_SIMGRID_LIBDIR}")
  elseif(ENV_SIMGRID_DIR)
    list(APPEND _lib_env "${ENV_SIMGRID_DIR}")
    list(APPEND _lib_env "${ENV_SIMGRID_DIR}/lib")
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

  # Try to find the simgrid lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(SIMGRID_LIBDIR)
    set(SIMGRID_simgrid_LIBRARY "SIMGRID_simgrid_LIBRARY-NOTFOUND")
    find_library(SIMGRID_simgrid_LIBRARY
      NAMES simgrid
      HINTS ${SIMGRID_LIBDIR})
  else()
    if(SIMGRID_DIR)
      set(SIMGRID_simgrid_LIBRARY "SIMGRID_simgrid_LIBRARY-NOTFOUND")
      find_library(SIMGRID_simgrid_LIBRARY
        NAMES simgrid
        HINTS ${SIMGRID_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(SIMGRID_simgrid_LIBRARY "SIMGRID_simgrid_LIBRARY-NOTFOUND")
      find_library(SIMGRID_simgrid_LIBRARY
        NAMES simgrid
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(SIMGRID_simgrid_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (SIMGRID_simgrid_LIBRARY)
    get_filename_component(simgrid_lib_path ${SIMGRID_simgrid_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(SIMGRID_LIBRARIES    "${SIMGRID_simgrid_LIBRARY}")
    set(SIMGRID_LIBRARY_DIRS "${simgrid_lib_path}")
  else ()
    set(SIMGRID_LIBRARIES    "SIMGRID_LIBRARIES-NOTFOUND")
    set(SIMGRID_LIBRARY_DIRS "SIMGRID_LIBRARY_DIRS-NOTFOUND")
    if(NOT SIMGRID_FIND_QUIETLY)
      message(STATUS "Looking for simgrid -- lib simgrid not found")
    endif()
  endif ()

  if (SIMGRID_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES SIMGRID_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT SIMGRID_FOUND) OR (SIMGRID_GIVEN_BY_USER) )

# check a function to validate the find
if(SIMGRID_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # SIMGRID
  if (SIMGRID_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${SIMGRID_INCLUDE_DIRS}")
  endif()
  if (SIMGRID_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${SIMGRID_CFLAGS_OTHER}")
  endif()
  if (SIMGRID_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${SIMGRID_LDFLAGS_OTHER}")
  endif()
  if (SIMGRID_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${SIMGRID_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${SIMGRID_LIBRARIES}")

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
  unset(SIMGRID_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(MSG_main SIMGRID_WORKS)
  mark_as_advanced(SIMGRID_WORKS)

  if(NOT SIMGRID_WORKS)
    if(NOT SIMGRID_FIND_QUIETLY)
      message(STATUS "Looking for simgrid : test of fut_keychange with simgrid library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(SIMGRID_LIBRARIES)

if (SIMGRID_LIBRARIES)
  if (SIMGRID_LIBRARY_DIRS)
    list(GET SIMGRID_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET SIMGRID_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(SIMGRID_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(SIMGRID_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of SIMGRID library" FORCE)
  else()
    set(SIMGRID_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of SIMGRID library" FORCE)
  endif()
endif()
mark_as_advanced(SIMGRID_DIR)
mark_as_advanced(SIMGRID_DIR_FOUND)

# check that SIMGRID has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SIMGRID DEFAULT_MSG
  SIMGRID_LIBRARIES
  SIMGRID_WORKS)
