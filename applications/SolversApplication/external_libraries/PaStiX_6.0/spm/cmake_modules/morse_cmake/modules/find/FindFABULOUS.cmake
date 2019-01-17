###
#
# @copyright (c) 2012-2018 Inria. All rights reserved.
#
###
#
# - Find FABULOUS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(FABULOUS
#               [REQUIRED]) # Fail with error if fabulous is not found
#
#  FABULOUS depends on the following libraries:
#   - CBLAS
#   - LAPACKE
#
# This module finds headers and fabulous library.
# Results are reported in variables:
#  FABULOUS_FOUND             - True if headers and requested libraries were found
#  FABULOUS_CFLAGS_OTHER      - fabulous compiler flags without headers paths
#  FABULOUS_LDFLAGS_OTHER     - fabulous linker flags without libraries
#  FABULOUS_INCLUDE_DIRS      - fabulous include directories
#  FABULOUS_MODULE_DIRS       - fabulous module directories for Fortran API
#  FABULOUS_LIBRARY_DIRS      - fabulous link directories
#  FABULOUS_LIBRARIES         - fabulous libraries to be linked (absolute path)
#  FABULOUS_CFLAGS_OTHER_DEP  - fabulous + dependencies compiler flags without headers paths
#  FABULOUS_LDFLAGS_OTHER_DEP - fabulous + dependencies linker flags without libraries
#  FABULOUS_INCLUDE_DIRS_DEP  - fabulous + dependencies include directories
#  FABULOUS_LIBRARY_DIRS_DEP  - fabulous + dependencies link directories
#  FABULOUS_LIBRARIES_DEP     - fabulous + dependencies libraries
#
#  FABULOUS_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = FABULOUS
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
# options at configure (ex: cmake path/to/project -DFABULOUS_DIR=path/to/fabulous):
#  FABULOUS_DIR             - Where to find the base directory of fabulous
#  FABULOUS_INCDIR          - Where to find the header files
#  FABULOUS_LIBDIR          - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: FABULOUS_DIR, FABULOUS_INCDIR, FABULOUS_LIBDIR

#=============================================================================
# Copyright 2012-2018 Inria
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

if (NOT FABULOUS_FOUND)
  set(FABULOUS_DIR "" CACHE PATH "Installation directory of FABULOUS library")
  if (NOT FABULOUS_FIND_QUIETLY)
    message(STATUS "A cache variable, namely FABULOUS_DIR, has been set to specify the install directory of FABULOUS")
  endif()
endif()

# FABULOUS depends on CBLAS
#---------------------------
if (NOT FABULOUS_FIND_QUIETLY)
  message(STATUS "Looking for FABULOUS - Try to detect CBLAS (depends on BLAS)")
endif()
if (FABULOUS_FIND_REQUIRED)
  find_package(CBLAS REQUIRED)
else()
  find_package(CBLAS)
endif()

# FABULOUS depends on LAPACKE
#-----------------------------
if (NOT FABULOUS_FIND_QUIETLY)
  message(STATUS "Looking for FABULOUS - Try to detect LAPACKE (depends on LAPACK)")
endif()
if (FABULOUS_FIND_REQUIRED)
  find_package(LAPACKE REQUIRED)
else()
  find_package(LAPACKE)
endif()
  
set(ENV_FABULOUS_DIR "$ENV{FABULOUS_DIR}")
set(ENV_FABULOUS_INCDIR "$ENV{FABULOUS_INCDIR}")
set(ENV_FABULOUS_LIBDIR "$ENV{FABULOUS_LIBDIR}")
set(FABULOUS_GIVEN_BY_USER "FALSE")
if ( FABULOUS_DIR OR ( FABULOUS_INCDIR AND FABULOUS_LIBDIR) OR ENV_FABULOUS_DIR OR (ENV_FABULOUS_INCDIR AND ENV_FABULOUS_LIBDIR) )
  set(FABULOUS_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if( PKG_CONFIG_EXECUTABLE AND NOT FABULOUS_GIVEN_BY_USER )

  pkg_search_module(FABULOUS fabulous_c_api)

  if (NOT FABULOUS_FIND_QUIETLY)
    if (FABULOUS_FOUND AND FABULOUS_LIBRARIES)
      message(STATUS "Looking for FABULOUS - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for FABULOUS - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing fabulous.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (FABULOUS_FOUND AND FABULOUS_LIBRARIES)
    set(FABULOUS_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(FABULOUS)
  else()
    set(FABULOUS_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif( PKG_CONFIG_EXECUTABLE AND NOT FABULOUS_GIVEN_BY_USER )

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT FABULOUS_FOUND) OR (FABULOUS_GIVEN_BY_USER) )

  if (NOT FABULOUS_FIND_QUIETLY)
    message(STATUS "Looking for FABULOUS - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  if(ENV_FABULOUS_INCDIR)
    list(APPEND _inc_env "${ENV_FABULOUS_INCDIR}")
  elseif(ENV_FABULOUS_DIR)
    list(APPEND _inc_env "${ENV_FABULOUS_DIR}")
    list(APPEND _inc_env "${ENV_FABULOUS_DIR}/include")
    list(APPEND _inc_env "${ENV_FABULOUS_DIR}/include/fabulous")
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

  # List of files to find
  set(FABULOUS_hdrs_to_find "fabulous.h;fabulous_mod.mod")

  # Try to find the fabulous header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(FABULOUS_INCDIR)
    foreach(_file ${FABULOUS_hdrs_to_find})
      set(FABULOUS_${_file}_DIRS "FABULOUS_${_file}_DIRS-NOTFOUND")
      find_path(FABULOUS_${_file}_DIRS
        NAMES ${_file}
        HINTS "${FABULOUS_INCDIR}" "${FABULOUS_INCDIR}/fabulous")
      mark_as_advanced(FABULOUS_${_file}_DIRS)
    endforeach()
  else()
    if(FABULOUS_DIR)
      foreach(_file ${FABULOUS_hdrs_to_find})
        set(FABULOUS_${_file}_DIRS "FABULOUS_${_file}_DIRS-NOTFOUND")
        find_path(FABULOUS_${_file}_DIRS
          NAMES ${_file}
          HINTS ${FABULOUS_DIR}
          PATH_SUFFIXES "include" "include/fabulous")
        mark_as_advanced(FABULOUS_${_file}_DIRS)
      endforeach()
    else()
      foreach(_file ${FABULOUS_hdrs_to_find})
        set(FABULOUS_${_file}_DIRS "FABULOUS_${_file}_DIRS-NOTFOUND")
        find_path(FABULOUS_${_file}_DIRS
          NAMES ${_file}
          HINTS ${PATH_TO_LOOK_FOR}
          PATH_SUFFIXES "fabulous")
        mark_as_advanced(FABULOUS_${_file}_DIRS)
      endforeach()
    endif()
  endif()

  # Add path to cmake variable
  # ------------------------------------
  if (FABULOUS_fabulous.h_DIRS)
    set(FABULOUS_INCLUDE_DIRS "${FABULOUS_fabulous.h_DIRS}")
  else ()
    set(FABULOUS_INCLUDE_DIRS "FABULOUS_INCLUDE_DIRS-NOTFOUND")
    if(NOT FABULOUS_FIND_QUIETLY)
      message(STATUS "Looking for fabulous -- fabulous.h not found")
    endif()
  endif ()
  if (FABULOUS_fabulous_mod.mod_DIRS)
    set(FABULOUS_MODULE_DIRS "${FABULOUS_fabulous_mod.mod_DIRS}")
  else ()
    set(FABULOUS_MODULE_DIRS "FABULOUS_MODULE_DIRS-NOTFOUND")
    if(NOT FABULOUS_FIND_QUIETLY)
      message(STATUS "Looking for fabulous -- fabulous_mod.mod not found")
    endif()
  endif ()

  if (FABULOUS_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES FABULOUS_INCLUDE_DIRS)
  endif ()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  if(ENV_FABULOUS_LIBDIR)
    list(APPEND _lib_env "${ENV_FABULOUS_LIBDIR}")
  elseif(ENV_FABULOUS_DIR)
    list(APPEND _lib_env "${ENV_FABULOUS_DIR}")
    list(APPEND _lib_env "${ENV_FABULOUS_DIR}/lib")
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

  # Try to find the fabulous lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(FABULOUS_LIBDIR)
    set(FABULOUS_fabulous_LIBRARY "FABULOUS_fabulous_LIBRARY-NOTFOUND")
    find_library(FABULOUS_fabulous_LIBRARY
      NAMES fabulous
      HINTS ${FABULOUS_LIBDIR})
  else()
    if(FABULOUS_DIR)
      set(FABULOUS_fabulous_LIBRARY "FABULOUS_fabulous_LIBRARY-NOTFOUND")
      find_library(FABULOUS_fabulous_LIBRARY
        NAMES fabulous
        HINTS ${FABULOUS_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(FABULOUS_fabulous_LIBRARY "FABULOUS_fabulous_LIBRARY-NOTFOUND")
      find_library(FABULOUS_fabulous_LIBRARY
        NAMES fabulous
        HINTS ${PATH_TO_LOOK_FOR})
    endif()
  endif()
  mark_as_advanced(FABULOUS_fabulous_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (FABULOUS_fabulous_LIBRARY)
    get_filename_component(fabulous_lib_path ${FABULOUS_fabulous_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    set(FABULOUS_LIBRARIES    "${FABULOUS_fabulous_LIBRARY}")
    set(FABULOUS_LIBRARY_DIRS "${fabulous_lib_path}")
  else ()
    set(FABULOUS_LIBRARIES    "FABULOUS_LIBRARIES-NOTFOUND")
    set(FABULOUS_LIBRARY_DIRS "FABULOUS_LIBRARY_DIRS-NOTFOUND")
    if(NOT FABULOUS_FIND_QUIETLY)
      message(STATUS "Looking for fabulous -- lib fabulous not found")
    endif()
  endif ()

  if (FABULOUS_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES FABULOUS_LIBRARY_DIRS)
  endif ()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT FABULOUS_FOUND) OR (FABULOUS_GIVEN_BY_USER) )

# check a function to validate the find
if(FABULOUS_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # FABULOUS
  if (FABULOUS_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${FABULOUS_INCLUDE_DIRS}")
  endif()
  if (FABULOUS_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${FABULOUS_CFLAGS_OTHER}")
  endif()
  if (FABULOUS_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${FABULOUS_LDFLAGS_OTHER}")
  endif()
  if (FABULOUS_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${FABULOUS_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${FABULOUS_LIBRARIES}")
  # LAPACKE
  if (LAPACKE_FOUND)
    if (LAPACKE_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${LAPACKE_INCLUDE_DIRS_DEP}")
    endif()
    if (LAPACKE_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${LAPACKE_CFLAGS_OTHER_DEP}")
    endif()
    if (LAPACKE_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${LAPACKE_LDFLAGS_OTHER_DEP}")
    endif()
    if(LAPACKE_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${LAPACKE_LIBRARY_DIRS_DEP}")
    endif()
    if (LAPACKE_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${LAPACKE_LIBRARIES_DEP}")
    endif()
  endif()
  # CBLAS
  if (CBLAS_FOUND)
    if (CBLAS_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS_DEP}")
    endif()
    if (CBLAS_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${CBLAS_CFLAGS_OTHER_DEP}")
    endif()
    if (CBLAS_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${CBLAS_LDFLAGS_OTHER_DEP}")
    endif()
    if(CBLAS_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS_DEP}")
    endif()
    if (CBLAS_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES_DEP}")
    endif()
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
  unset(FABULOUS_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(fabulous_solve FABULOUS_WORKS)
  mark_as_advanced(FABULOUS_WORKS)

  if(FABULOUS_WORKS)
    # save link with dependencies
    set(FABULOUS_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(FABULOUS_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(FABULOUS_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(FABULOUS_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(FABULOUS_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT FABULOUS_FIND_QUIETLY)
      message(STATUS "Looking for fabulous : test of fabulous_topology_init with fabulous library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(FABULOUS_LIBRARIES)

if (FABULOUS_LIBRARIES)
  if (FABULOUS_LIBRARY_DIRS)
    list(GET FABULOUS_LIBRARY_DIRS 0 first_lib_path)
  else()
    list(GET FABULOUS_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(FABULOUS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(FABULOUS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of FABULOUS library" FORCE)
  else()
    set(FABULOUS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of FABULOUS library" FORCE)
  endif()
endif()
mark_as_advanced(FABULOUS_DIR)
mark_as_advanced(FABULOUS_DIR_FOUND)

# check that FABULOUS has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FABULOUS DEFAULT_MSG
  FABULOUS_LIBRARIES
  FABULOUS_WORKS)

