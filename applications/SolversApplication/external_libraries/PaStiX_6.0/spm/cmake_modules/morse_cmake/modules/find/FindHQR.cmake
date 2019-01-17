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
# - Find HQR include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(HQR
#               [REQUIRED]) # Fail with error if hqr is not found
#
# This module finds headers and hqr library.
# Results are reported in variables:
#  HQR_FOUND           - True if headers and requested libraries were found
#  HQR_INCLUDE_DIRS    - hqr include directories
#  HQR_LIBRARY_DIRS    - Link directories for hqr libraries
#  HQR_LIBRARIES       - hqr component libraries to be linked
#
#  HQR_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = HQR
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
# options at configure (ex: cmake path/to/project -DHQR_DIR=path/to/hqr):
#  HQR_DIR             - Where to find the base directory of hqr
#  HQR_INCDIR          - Where to find the header files
#  HQR_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: HQR_DIR, HQR_INCDIR, HQR_LIBDIR

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

if (NOT HQR_FOUND)
  set(HQR_DIR "" CACHE PATH "Installation directory of HQR library")
  if (NOT HQR_FIND_QUIETLY)
    message(STATUS "A cache variable, namely HQR_DIR, has been set to specify the install directory of HQR")
  endif()
endif()

set(ENV_HQR_DIR "$ENV{HQR_DIR}")
set(ENV_HQR_INCDIR "$ENV{HQR_INCDIR}")
set(ENV_HQR_LIBDIR "$ENV{HQR_LIBDIR}")
set(HQR_GIVEN_BY_USER "FALSE")
if ( HQR_DIR OR ( HQR_INCDIR AND HQR_LIBDIR) OR ENV_HQR_DIR OR (ENV_HQR_INCDIR AND ENV_HQR_LIBDIR) )
  set(HQR_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT HQR_GIVEN_BY_USER)

  pkg_search_module(HQR hqr)

  if (NOT HQR_FIND_QUIETLY)
    if (HQR_FOUND AND HQR_LIBRARIES)
      message(STATUS "Looking for HQR - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for HQR - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing hqr.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()
  if (HQR_FOUND AND HQR_LIBRARIES)
    set(HQR_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(HQR)
  else()
    set(HQR_FOUND_WITH_PKGCONFIG "FALSE")
  endif()
endif(PKG_CONFIG_EXECUTABLE AND NOT HQR_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT HQR_FOUND) OR (HQR_GIVEN_BY_USER) )

  if (NOT HQR_FIND_QUIETLY)
    message(STATUS "Looking for HQR - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_HQR_DIR "$ENV{HQR_DIR}")
  set(ENV_HQR_INCDIR "$ENV{HQR_INCDIR}")
  if(ENV_HQR_INCDIR)
    list(APPEND _inc_env "${ENV_HQR_INCDIR}")
  elseif(ENV_HQR_DIR)
    list(APPEND _inc_env "${ENV_HQR_DIR}")
    list(APPEND _inc_env "${ENV_HQR_DIR}/include")
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

  # Try to find the hqr header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(HQR_INCDIR)
    set(HQR_HQR.h_DIRS "HQR_libhqr.h_DIRS-NOTFOUND")
    find_path(HQR_libhqr.h_DIRS
      NAMES libhqr.h
      HINTS ${HQR_INCDIR})
  else()
    if(HQR_DIR)
      set(HQR_libhqr.h_DIRS "HQR_libhqr.h_DIRS-NOTFOUND")
      find_path(HQR_libhqr.h_DIRS
        NAMES libhqr.h
        HINTS ${HQR_DIR})
    else()
      set(HQR_libhqr.h_DIRS "HQR_libhqr.h_DIRS-NOTFOUND")
      find_path(HQR_libhqr.h_DIRS
        NAMES libhqr.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "hqr")
    endif()
  endif()
  mark_as_advanced(HQR_libhqr.h_DIRS)

  # Add path to cmake variable
  # ------------------------------------
  if (HQR_libhqr.h_DIRS)
    set(HQR_INCLUDE_DIRS "${HQR_libhqr.h_DIRS}")
  else ()
    set(HQR_INCLUDE_DIRS "HQR_INCLUDE_DIRS-NOTFOUND")
    if(NOT HQR_FIND_QUIETLY)
      message(STATUS "Looking for hqr -- libhqr.h not found")
    endif()
  endif ()

  if (HQR_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES HQR_INCLUDE_DIRS)
  endif ()

  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_HQR_LIBDIR "$ENV{HQR_LIBDIR}")
  if(ENV_HQR_LIBDIR)
    list(APPEND _lib_env "${ENV_HQR_LIBDIR}")
  elseif(ENV_HQR_DIR)
    list(APPEND _lib_env "${ENV_HQR_DIR}")
    list(APPEND _lib_env "${ENV_HQR_DIR}/lib")
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

  # Try to find the hqr lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(HQR_LIBDIR)
    set(HQR_hqr_LIBRARY "HQR_hqr_LIBRARY-NOTFOUND")
    find_library(HQR_hqr_LIBRARY
      NAMES hqr
      HINTS ${HQR_LIBDIR})
  else()
    if(HQR_DIR)
      set(HQR_hqr_LIBRARY "HQR_hqr_LIBRARY-NOTFOUND")
      find_library(HQR_hqr_LIBRARY
        NAMES hqr
        HINTS ${HQR_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(HQR_hqr_LIBRARY "HQR_hqr_LIBRARY-NOTFOUND")
      find_library(HQR_hqr_LIBRARY
        NAMES hqr
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(HQR_hqr_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (HQR_hqr_LIBRARY)
    get_filename_component(hqr_lib_path ${HQR_hqr_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(HQR_LIBRARIES    "${HQR_hqr_LIBRARY}")
    set(HQR_LIBRARY_DIRS "${hqr_lib_path}")
  else ()
    set(HQR_LIBRARIES    "HQR_LIBRARIES-NOTFOUND")
    set(HQR_LIBRARY_DIRS "HQR_LIBRARY_DIRS-NOTFOUND")
    if(NOT HQR_FIND_QUIETLY)
      message(STATUS "Looking for hqr -- lib hqr not found")
    endif()
  endif ()

  if (HQR_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES HQR_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT HQR_FOUND) OR (HQR_GIVEN_BY_USER) )

# check a function to validate the find
if(HQR_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # HQR
  if (HQR_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${HQR_INCLUDE_DIRS}")
  endif()
  if (HQR_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${HQR_CFLAGS_OTHER}")
  endif()
  if (HQR_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${HQR_LDFLAGS_OTHER}")
  endif()
  if (HQR_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${HQR_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${HQR_LIBRARIES}")

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
  unset(HQR_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(libhqr_init_hqr HQR_WORKS)
  mark_as_advanced(HQR_WORKS)

  if(NOT HQR_WORKS)
    if(NOT HQR_FIND_QUIETLY)
      message(STATUS "Looking for hqr : test of libhqr_hqr_init with hqr library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(HQR_LIBRARIES)

if (HQR_LIBRARIES)
  if (HQR_LIBRARY_DIRS)
    list(GET HQR_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET HQR_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(HQR_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(HQR_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of HQR library" FORCE)
  else()
    set(HQR_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of HQR library" FORCE)
  endif()
endif()
mark_as_advanced(HQR_DIR)
mark_as_advanced(HQR_DIR_FOUND)

# check that HQR has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(HQR DEFAULT_MSG
  HQR_LIBRARIES
  HQR_WORKS)
