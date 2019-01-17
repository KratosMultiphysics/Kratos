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
# - Find EZTRACE include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(EZTRACE
#               [REQUIRED]) # Fail with error if eztrace is not found
#
#  EZTRACE depends on the following libraries:
#   - libbfd.so
#   - libiberty.so (optional)
#
# This module finds headers and eztrace library.
# Results are reported in variables:
#  EZTRACE_FOUND           - True if headers and requested libraries were found
#  EZTRACE_INCLUDE_DIRS    - eztrace include directories
#  EZTRACE_LIBRARY_DIRS    - Link directories for eztrace libraries
#  EZTRACE_LIBRARIES       - eztrace component libraries to be linked
#
#  EZTRACE_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = EZTRACE
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
# options at configure (ex: cmake path/to/project -DEZTRACE_DIR=path/to/eztrace):
#  EZTRACE_DIR             - Where to find the base directory of eztrace
#  EZTRACE_INCDIR          - Where to find the header files
#  EZTRACE_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: EZTRACE_DIR, EZTRACE_INCDIR, EZTRACE_LIBDIR

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

if (NOT EZTRACE_FOUND)
  set(EZTRACE_DIR "" CACHE PATH "Installation directory of EZTRACE library")
  if (NOT EZTRACE_FIND_QUIETLY)
    message(STATUS "A cache variable, namely EZTRACE_DIR, has been set to specify the install directory of EZTRACE")
  endif()
endif()

set(ENV_EZTRACE_DIR "$ENV{EZTRACE_DIR}")
set(ENV_EZTRACE_INCDIR "$ENV{EZTRACE_INCDIR}")
set(ENV_EZTRACE_LIBDIR "$ENV{EZTRACE_LIBDIR}")
set(EZTRACE_GIVEN_BY_USER "FALSE")
if ( EZTRACE_DIR OR ( EZTRACE_INCDIR AND EZTRACE_LIBDIR) OR ENV_EZTRACE_DIR OR (ENV_EZTRACE_INCDIR AND ENV_EZTRACE_LIBDIR) )
  set(EZTRACE_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if( PKG_CONFIG_EXECUTABLE AND NOT EZTRACE_GIVEN_BY_USER )

  pkg_search_module(EZTRACE eztrace)

  if (NOT EZTRACE_FIND_QUIETLY)
    if (EZTRACE_FOUND AND EZTRACE_LIBRARIES)
      message(STATUS "Looking for EZTRACE - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for EZTRACE - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing eztrace.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
    if (EZTRACE_FOUND AND EZTRACE_LIBRARIES)
      set(EZTRACE_FOUND_WITH_PKGCONFIG "TRUE")
      find_pkgconfig_libraries_absolute_path(EZTRACE)
    else()
      set(EZTRACE_FOUND_WITH_PKGCONFIG "FALSE")
    endif()
  endif()

endif( PKG_CONFIG_EXECUTABLE AND NOT EZTRACE_GIVEN_BY_USER )

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT EZTRACE_FOUND) OR (EZTRACE_GIVEN_BY_USER) )

  if (NOT EZTRACE_FIND_QUIETLY)
    message(STATUS "Looking for EZTRACE - PkgConfig not used")
  endif()

  # Looking for libbfd
  # ------------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
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
  list(REMOVE_DUPLICATES _lib_env)

  # set paths where to look for
  set(PATH_TO_LOOK_FOR "${_lib_env}")

  # Try to find the bfd lib in the given paths
  # ------------------------------------------

  # call cmake macro to find the lib path
  set(EZTRACE_bfd_LIBRARY "EZTRACE_bfd_LIBRARY-NOTFOUND")
  find_library(EZTRACE_bfd_LIBRARY
    NAMES bfd
    HINTS ${PATH_TO_LOOK_FOR})
  mark_as_advanced(EZTRACE_bfd_LIBRARY)

  if (NOT EZTRACE_bfd_LIBRARY)
    if(NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "Looking for eztrace -- lib bfd not found")
    endif()
  endif ()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  if(ENV_EZTRACE_INCDIR)
    list(APPEND _inc_env "${ENV_EZTRACE_INCDIR}")
  elseif(ENV_EZTRACE_DIR)
    list(APPEND _inc_env "${ENV_EZTRACE_DIR}")
    list(APPEND _inc_env "${ENV_EZTRACE_DIR}/include")
    list(APPEND _inc_env "${ENV_EZTRACE_DIR}/include/eztrace")
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

  # Try to find the eztrace header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(EZTRACE_INCDIR)
    set(EZTRACE_eztrace.h_DIRS "EZTRACE_eztrace.h_DIRS-NOTFOUND")
    find_path(EZTRACE_eztrace.h_DIRS
      NAMES eztrace.h
      HINTS ${EZTRACE_INCDIR})
  else()
    if(EZTRACE_DIR)
      set(EZTRACE_eztrace.h_DIRS "EZTRACE_eztrace.h_DIRS-NOTFOUND")
      find_path(EZTRACE_eztrace.h_DIRS
        NAMES eztrace.h
        HINTS ${EZTRACE_DIR}
        PATH_SUFFIXES "include" "include/eztrace")
    else()
      set(EZTRACE_eztrace.h_DIRS "EZTRACE_eztrace.h_DIRS-NOTFOUND")
      find_path(EZTRACE_eztrace.h_DIRS
        NAMES eztrace.h
        HINTS ${PATH_TO_LOOK_FOR}
        PATH_SUFFIXES "eztrace")
    endif()
  endif()
  mark_as_advanced(EZTRACE_eztrace.h_DIRS)

  # Add path to cmake variable
  # ------------------------------------
  if (EZTRACE_eztrace.h_DIRS)
    set(EZTRACE_INCLUDE_DIRS "${EZTRACE_eztrace.h_DIRS}")
  else ()
    set(EZTRACE_INCLUDE_DIRS "EZTRACE_INCLUDE_DIRS-NOTFOUND")
    if(NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "Looking for eztrace -- eztrace.h not found")
    endif()
  endif ()

  if (EZTRACE_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES EZTRACE_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  if(ENV_EZTRACE_LIBDIR)
    list(APPEND _lib_env "${ENV_EZTRACE_LIBDIR}")
  elseif(ENV_EZTRACE_DIR)
    list(APPEND _lib_env "${ENV_EZTRACE_DIR}")
    list(APPEND _lib_env "${ENV_EZTRACE_DIR}/lib")
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

  # Try to find the eztrace lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(EZTRACE_LIBDIR)
    set(EZTRACE_eztrace_LIBRARY "EZTRACE_eztrace_LIBRARY-NOTFOUND")
    find_library(EZTRACE_eztrace_LIBRARY
      NAMES eztrace
      HINTS ${EZTRACE_LIBDIR})
  else()
    if(EZTRACE_DIR)
      set(EZTRACE_eztrace_LIBRARY "EZTRACE_eztrace_LIBRARY-NOTFOUND")
      find_library(EZTRACE_eztrace_LIBRARY
        NAMES eztrace
        HINTS ${EZTRACE_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(EZTRACE_eztrace_LIBRARY "EZTRACE_eztrace_LIBRARY-NOTFOUND")
      find_library(EZTRACE_eztrace_LIBRARY
        NAMES eztrace
        HINTS ${PATH_TO_LOOK_FOR})
    endif()
  endif()
  mark_as_advanced(EZTRACE_eztrace_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (EZTRACE_eztrace_LIBRARY)
    get_filename_component(eztrace_lib_path ${EZTRACE_eztrace_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(EZTRACE_LIBRARIES    "${EZTRACE_eztrace_LIBRARY}")
    set(EZTRACE_LIBRARY_DIRS "${eztrace_lib_path}")
  else ()
    set(EZTRACE_LIBRARIES    "EZTRACE_LIBRARIES-NOTFOUND")
    set(EZTRACE_LIBRARY_DIRS "EZTRACE_LIBRARY_DIRS-NOTFOUND")
    if(NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "Looking for eztrace -- lib eztrace not found")
    endif()
  endif ()

  if (EZTRACE_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES EZTRACE_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT EZTRACE_FOUND) OR (EZTRACE_GIVEN_BY_USER) )

# check a function to validate the find
if(EZTRACE_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # EZTRACE
  if (EZTRACE_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${EZTRACE_INCLUDE_DIRS}")
  endif()
  if (EZTRACE_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${EZTRACE_CFLAGS_OTHER}")
  endif()
  if (EZTRACE_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${EZTRACE_LDFLAGS_OTHER}")
  endif()
  if (EZTRACE_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${EZTRACE_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${EZTRACE_LIBRARIES}")

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
  unset(EZTRACE_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(eztrace_start EZTRACE_WORKS)
  mark_as_advanced(EZTRACE_WORKS)

  if(NOT EZTRACE_WORKS)
    if(NOT EZTRACE_FIND_QUIETLY)
      message(STATUS "Looking for eztrace : test of eztrace_topology_init with eztrace library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(EZTRACE_LIBRARIES)

if (EZTRACE_LIBRARIES)
  if (EZTRACE_LIBRARY_DIRS)
    list(GET EZTRACE_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET EZTRACE_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(EZTRACE_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(EZTRACE_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of EZTRACE library" FORCE)
  else()
    set(EZTRACE_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of EZTRACE library" FORCE)
  endif()
endif()
mark_as_advanced(EZTRACE_DIR)
mark_as_advanced(EZTRACE_DIR_FOUND)

# check that EZTRACE has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(EZTRACE DEFAULT_MSG
  EZTRACE_LIBRARIES
  EZTRACE_WORKS)
