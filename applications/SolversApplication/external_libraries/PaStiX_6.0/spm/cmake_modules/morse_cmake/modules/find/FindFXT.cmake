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
# - Find FXT include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(FXT
#               [REQUIRED]) # Fail with error if fxt is not found
#
# This module finds headers and fxt library.
# Results are reported in variables:
#  FXT_FOUND           - True if headers and requested libraries were found
#  FXT_INCLUDE_DIRS    - fxt include directories
#  FXT_LIBRARY_DIRS    - fxt link directories
#  FXT_LIBRARIES       - fxt libraries to be linked (absolute path)
#
#  FXT_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = FXT
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
# options at configure (ex: cmake path/to/project -DFXT_DIR=path/to/fxt):
#  FXT_DIR             - Where to find the base directory of fxt
#  FXT_INCDIR          - Where to find the header files
#  FXT_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: FXT_DIR, FXT_INCDIR, FXT_LIBDIR

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

if (NOT FXT_FOUND)
  set(FXT_DIR "" CACHE PATH "Installation directory of FXT library")
  if (NOT FXT_FIND_QUIETLY)
    message(STATUS "A cache variable, namely FXT_DIR, has been set to specify the install directory of FXT")
  endif()
endif()

set(ENV_FXT_DIR "$ENV{FXT_DIR}")
set(ENV_FXT_INCDIR "$ENV{FXT_INCDIR}")
set(ENV_FXT_LIBDIR "$ENV{FXT_LIBDIR}")
set(FXT_GIVEN_BY_USER "FALSE")
if ( FXT_DIR OR ( FXT_INCDIR AND FXT_LIBDIR) OR ENV_FXT_DIR OR (ENV_FXT_INCDIR AND ENV_FXT_LIBDIR) )
  set(FXT_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT FXT_GIVEN_BY_USER)

  pkg_search_module(FXT fxt)

  if (NOT FXT_FIND_QUIETLY)
    if (FXT_FOUND AND FXT_LIBRARIES)
      message(STATUS "Looking for FXT - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for FXT - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing fxt.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (FXT_FOUND AND FXT_LIBRARIES)
    set(FXT_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(FXT)
  else()
    set(FXT_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT FXT_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT FXT_FOUND) OR (FXT_GIVEN_BY_USER) )

  if (NOT FXT_FIND_QUIETLY)
    message(STATUS "Looking for FXT - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_FXT_DIR "$ENV{FXT_DIR}")
  set(ENV_FXT_INCDIR "$ENV{FXT_INCDIR}")
  if(ENV_FXT_INCDIR)
    list(APPEND _inc_env "${ENV_FXT_INCDIR}")
  elseif(ENV_FXT_DIR)
    list(APPEND _inc_env "${ENV_FXT_DIR}")
    list(APPEND _inc_env "${ENV_FXT_DIR}/include")
    list(APPEND _inc_env "${ENV_FXT_DIR}/include/fxt")
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

  # Try to find the fxt header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(FXT_INCDIR)
    set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
    find_path(FXT_fxt.h_DIRS
      NAMES fxt.h
      HINTS ${FXT_INCDIR})
  else()
    if(FXT_DIR)
      set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
      find_path(FXT_fxt.h_DIRS
        NAMES fxt.h
        HINTS ${FXT_DIR}
        PATH_SUFFIXES "include" "include/fxt")
    else()
      set(FXT_fxt.h_DIRS "FXT_fxt.h_DIRS-NOTFOUND")
      find_path(FXT_fxt.h_DIRS
        NAMES fxt.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "fxt")
    endif()
  endif()
  mark_as_advanced(FXT_fxt.h_DIRS)

  # Add path to cmake variable
  # ------------------------------------
  if (FXT_fxt.h_DIRS)
    set(FXT_INCLUDE_DIRS "${FXT_fxt.h_DIRS}")
  else ()
    set(FXT_INCLUDE_DIRS "FXT_INCLUDE_DIRS-NOTFOUND")
    if(NOT FXT_FIND_QUIETLY)
      message(STATUS "Looking for fxt -- fxt.h not found")
    endif()
  endif ()

  if (FXT_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES FXT_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_FXT_LIBDIR "$ENV{FXT_LIBDIR}")
  if(ENV_FXT_LIBDIR)
    list(APPEND _lib_env "${ENV_FXT_LIBDIR}")
  elseif(ENV_FXT_DIR)
    list(APPEND _lib_env "${ENV_FXT_DIR}")
    list(APPEND _lib_env "${ENV_FXT_DIR}/lib")
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

  # Try to find the fxt lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(FXT_LIBDIR)
    set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
    find_library(FXT_fxt_LIBRARY
      NAMES fxt
      HINTS ${FXT_LIBDIR})
  else()
    if(FXT_DIR)
      set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
      find_library(FXT_fxt_LIBRARY
        NAMES fxt
        HINTS ${FXT_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(FXT_fxt_LIBRARY "FXT_fxt_LIBRARY-NOTFOUND")
      find_library(FXT_fxt_LIBRARY
        NAMES fxt
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(FXT_fxt_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (FXT_fxt_LIBRARY)
    get_filename_component(fxt_lib_path ${FXT_fxt_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(FXT_LIBRARIES    "${FXT_fxt_LIBRARY}")
    set(FXT_LIBRARY_DIRS "${fxt_lib_path}")
  else ()
    set(FXT_LIBRARIES    "FXT_LIBRARIES-NOTFOUND")
    set(FXT_LIBRARY_DIRS "FXT_LIBRARY_DIRS-NOTFOUND")
    if(NOT FXT_FIND_QUIETLY)
      message(STATUS "Looking for fxt -- lib fxt not found")
    endif()
  endif ()

  if (FXT_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES FXT_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT FXT_FOUND) OR (FXT_GIVEN_BY_USER) )

# check a function to validate the find
if(FXT_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # FXT
  if (FXT_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${FXT_INCLUDE_DIRS}")
  endif()
  if (FXT_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${FXT_CFLAGS_OTHER}")
  endif()
  if (FXT_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${FXT_LDFLAGS_OTHER}")
  endif()
  if (FXT_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${FXT_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${FXT_LIBRARIES}")

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
  unset(FXT_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(fut_keychange FXT_WORKS)
  mark_as_advanced(FXT_WORKS)

  if(NOT FXT_WORKS)
    if(NOT FXT_FIND_QUIETLY)
      message(STATUS "Looking for fxt : test of fut_keychange with fxt library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(FXT_LIBRARIES)

if (FXT_LIBRARIES)
  if (FXT_LIBRARY_DIRS)
    list(GET FXT_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET FXT_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(FXT_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(FXT_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of FXT library" FORCE)
  else()
    set(FXT_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of FXT library" FORCE)
  endif()
endif()
mark_as_advanced(FXT_DIR)
mark_as_advanced(FXT_DIR_FOUND)

# check that FXT has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FXT DEFAULT_MSG
  FXT_LIBRARIES
  FXT_WORKS)
