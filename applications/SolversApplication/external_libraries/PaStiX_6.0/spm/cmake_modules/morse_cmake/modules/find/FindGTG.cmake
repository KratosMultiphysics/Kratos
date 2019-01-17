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
# - Find GTG include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(GTG
#               [REQUIRED]) # Fail with error if gtg is not found
#
# This module finds headers and gtg library.
# Results are reported in variables:
#  GTG_FOUND           - True if headers and requested libraries were found
#  GTG_INCLUDE_DIRS    - gtg include directories
#  GTG_LIBRARY_DIRS    - Link directories for gtg libraries
#  GTG_LIBRARIES       - gtg component libraries to be linked
#
#  GTG_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = GTG
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
# options at configure (ex: cmake path/to/project -DGTG_DIR=path/to/gtg):
#  GTG_DIR             - Where to find the base directory of gtg
#  GTG_INCDIR          - Where to find the header files
#  GTG_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: GTG_DIR, GTG_INCDIR, GTG_LIBDIR

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

if (NOT GTG_FOUND)
  set(GTG_DIR "" CACHE PATH "Installation directory of GTG library")
  if (NOT GTG_FIND_QUIETLY)
    message(STATUS "A cache variable, namely GTG_DIR, has been set to specify the install directory of GTG")
  endif()
endif()

set(ENV_GTG_DIR "$ENV{GTG_DIR}")
set(ENV_GTG_INCDIR "$ENV{GTG_INCDIR}")
set(ENV_GTG_LIBDIR "$ENV{GTG_LIBDIR}")
set(GTG_GIVEN_BY_USER "FALSE")
if ( GTG_DIR OR ( GTG_INCDIR AND GTG_LIBDIR) OR ENV_GTG_DIR OR (ENV_GTG_INCDIR AND ENV_GTG_LIBDIR) )
  set(GTG_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT GTG_GIVEN_BY_USER)

  pkg_search_module(GTG gtg)

  if (NOT GTG_FIND_QUIETLY)
    if (GTG_FOUND AND GTG_LIBRARIES)
      message(STATUS "Looking for GTG - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for GTG - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing gtg.pc to the"
        "\n   PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (GTG_FOUND AND GTG_LIBRARIES)
    set(GTG_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(GTG)
  else()
    set(GTG_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT GTG_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT GTG_FOUND) OR (GTG_GIVEN_BY_USER) )

  if (NOT GTG_FIND_QUIETLY)
    message(STATUS "Looking for GTG - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_GTG_DIR "$ENV{GTG_DIR}")
  set(ENV_GTG_INCDIR "$ENV{GTG_INCDIR}")
  if(ENV_GTG_INCDIR)
    list(APPEND _inc_env "${ENV_GTG_INCDIR}")
  elseif(ENV_GTG_DIR)
    list(APPEND _inc_env "${ENV_GTG_DIR}")
    list(APPEND _inc_env "${ENV_GTG_DIR}/include")
    list(APPEND _inc_env "${ENV_GTG_DIR}/include/gtg")
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

  # Try to find the gtg header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(GTG_INCDIR)
    set(GTG_GTG.h_DIRS "GTG_GTG.h_DIRS-NOTFOUND")
    find_path(GTG_GTG.h_DIRS
      NAMES GTG.h
      HINTS ${GTG_INCDIR})
  else()
    if(GTG_DIR)
      set(GTG_GTG.h_DIRS "GTG_GTG.h_DIRS-NOTFOUND")
      find_path(GTG_GTG.h_DIRS
        NAMES GTG.h
        HINTS ${GTG_DIR}
        PATH_SUFFIXES "include" "include/gtg")
    else()
      set(GTG_GTG.h_DIRS "GTG_GTG.h_DIRS-NOTFOUND")
      find_path(GTG_GTG.h_DIRS
        NAMES GTG.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "gtg")
    endif()
  endif()
  mark_as_advanced(GTG_GTG.h_DIRS)

  # Add path to cmake variable
  # ------------------------------------
  if (GTG_GTG.h_DIRS)
    set(GTG_INCLUDE_DIRS "${GTG_GTG.h_DIRS}")
  else ()
    set(GTG_INCLUDE_DIRS "GTG_INCLUDE_DIRS-NOTFOUND")
    if(NOT GTG_FIND_QUIETLY)
      message(STATUS "Looking for gtg -- GTG.h not found")
    endif()
  endif ()

  if (GTG_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES GTG_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_GTG_LIBDIR "$ENV{GTG_LIBDIR}")
  if(ENV_GTG_LIBDIR)
    list(APPEND _lib_env "${ENV_GTG_LIBDIR}")
  elseif(ENV_GTG_DIR)
    list(APPEND _lib_env "${ENV_GTG_DIR}")
    list(APPEND _lib_env "${ENV_GTG_DIR}/lib")
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

  # Try to find the gtg lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(GTG_LIBDIR)
    set(GTG_gtg_LIBRARY "GTG_gtg_LIBRARY-NOTFOUND")
    find_library(GTG_gtg_LIBRARY
      NAMES gtg
      HINTS ${GTG_LIBDIR})
  else()
    if(GTG_DIR)
      set(GTG_gtg_LIBRARY "GTG_gtg_LIBRARY-NOTFOUND")
      find_library(GTG_gtg_LIBRARY
        NAMES gtg
        HINTS ${GTG_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(GTG_gtg_LIBRARY "GTG_gtg_LIBRARY-NOTFOUND")
      find_library(GTG_gtg_LIBRARY
        NAMES gtg
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(GTG_gtg_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (GTG_gtg_LIBRARY)
    get_filename_component(gtg_lib_path ${GTG_gtg_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(GTG_LIBRARIES    "${GTG_gtg_LIBRARY}")
    set(GTG_LIBRARY_DIRS "${gtg_lib_path}")
  else ()
    set(GTG_LIBRARIES    "GTG_LIBRARIES-NOTFOUND")
    set(GTG_LIBRARY_DIRS "GTG_LIBRARY_DIRS-NOTFOUND")
    if(NOT GTG_FIND_QUIETLY)
      message(STATUS "Looking for gtg -- lib gtg not found")
    endif()
  endif ()

  if (GTG_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES GTG_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT GTG_FOUND) OR (GTG_GIVEN_BY_USER) )

# check a function to validate the find
if(GTG_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # GTG
  if (GTG_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${GTG_INCLUDE_DIRS}")
  endif()
  if (GTG_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${GTG_CFLAGS_OTHER}")
  endif()
  if (GTG_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${GTG_LDFLAGS_OTHER}")
  endif()
  if (GTG_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${GTG_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${GTG_LIBRARIES}")

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
  unset(GTG_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(initTrace GTG_WORKS)
  mark_as_advanced(GTG_WORKS)

  if(NOT GTG_WORKS)
    if(NOT GTG_FIND_QUIETLY)
      message(STATUS "Looking for gtg : test of GTG_start with gtg library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(GTG_LIBRARIES)

if (GTG_LIBRARIES)
  if (GTG_LIBRARY_DIRS)
    list(GET GTG_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET GTG_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(GTG_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(GTG_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of GTG library" FORCE)
  else()
    set(GTG_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of GTG library" FORCE)
  endif()
endif()
mark_as_advanced(GTG_DIR)
mark_as_advanced(GTG_DIR_FOUND)

# check that GTG has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GTG DEFAULT_MSG
  GTG_LIBRARIES
  GTG_WORKS)
