###
#
# @copyright (c) 2018 Inria. All rights reserved.
#
###
#
# - Find SPM include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(SPM
#               [REQUIRED] # Fail with error if spm is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  SPM depends on the following libraries:
#   - CBLAS
#   - LAPACKE
#
#  COMPONENTS are options.
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - FORTRAN: to provide Fortran libraries of SPM in LIBRARIES
#
##
# This module finds headers and spm library.
# Results are reported in variables:
#  SPM_FOUND             - True if headers and requested libraries were found
#  SPM_CFLAGS_OTHER      - spm compiler flags without headers paths
#  SPM_LDFLAGS_OTHER     - spm linker flags without libraries
#  SPM_INCLUDE_DIRS      - spm include directories
#  SPM_LIBRARY_DIRS      - spm link directories
#  SPM_LIBRARIES         - spm libraries to be linked (absolute path)
#  SPM_LIBRARIES_FORTRAN - spm Fortran libraries to be linked (absolute path)
#  SPM_CFLAGS_OTHER_DEP  - spm + dependencies compiler flags without headers paths
#  SPM_LDFLAGS_OTHER_DEP - spm + dependencies linker flags without libraries
#  SPM_INCLUDE_DIRS_DEP  - spm + dependencies include directories
#  SPM_LIBRARY_DIRS_DEP  - spm + dependencies link directories
#  SPM_LIBRARIES_DEP     - spm + dependencies libraries
#  SPM_INTSIZE           - Number of octets occupied by a spm_int_t
#
#  SPM_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = SPM
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
# options at configure (ex: cmake path/to/project -DSPM_DIR=path/to/spm):
#  SPM_DIR              - Where to find the base directory of spm
#  SPM_INCDIR           - Where to find the header files
#  SPM_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: SPM_DIR, SPM_INCDIR, SPM_LIBDIR

#=============================================================================
# Copyright 2018 Inria
# Copyright 2018 Florent Pruvost
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

if (NOT SPM_FOUND)
  set(SPM_DIR "" CACHE PATH "Installation directory of SPM library")
  if (NOT SPM_FIND_QUIETLY)
    message(STATUS "A cache variable, namely SPM_DIR, has been set to specify the install directory of SPM")
  endif()
endif()

set(SPM_LOOK_FOR_FORTRAN OFF)
if( SPM_FIND_COMPONENTS )
  foreach( component ${SPM_FIND_COMPONENTS} )
    if (${component} STREQUAL "FORTRAN")
      set(SPM_LOOK_FOR_FORTRAN ON)
    endif()
  endforeach()
endif()

# Required dependencies
# ---------------------

# SPM depends on CBLAS
#-----------------------
if (NOT SPM_FIND_QUIETLY)
  message(STATUS "Looking for SPM - Try to detect CBLAS")
endif()
if (SPM_FIND_REQUIRED)
  find_package(CBLAS REQUIRED)
else()
  find_package(CBLAS)
endif()

# SPM depends on LAPACKE
#-----------------------
if (NOT SPM_FIND_QUIETLY)
  message(STATUS "Looking for SPM - Try to detect LAPACKE")
endif()
if (SPM_FIND_REQUIRED)
  find_package(LAPACKE REQUIRED)
else()
  find_package(LAPACKE)
endif()

set(ENV_SPM_DIR "$ENV{SPM_DIR}")
set(ENV_SPM_INCDIR "$ENV{SPM_INCDIR}")
set(ENV_SPM_LIBDIR "$ENV{SPM_LIBDIR}")
set(SPM_GIVEN_BY_USER "FALSE")
if ( SPM_DIR OR ( SPM_INCDIR AND SPM_LIBDIR) OR ENV_SPM_DIR OR (ENV_SPM_INCDIR AND ENV_SPM_LIBDIR) )
  set(SPM_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT SPM_GIVEN_BY_USER)

  pkg_search_module(SPM spm)

  if (NOT SPM_FIND_QUIETLY)
    if (SPM_FOUND AND SPM_LIBRARIES)
      message(STATUS "Looking for SPM - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for SPM - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing spm.pc"
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (SPM_FOUND AND SPM_LIBRARIES)
    set(SPM_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(SPM)
  else()
    set(SPM_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT SPM_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT SPM_FOUND) OR (SPM_GIVEN_BY_USER) )

  # Looking for SPM
  # ------------------

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_SPM_DIR "$ENV{SPM_DIR}")
  set(ENV_SPM_INCDIR "$ENV{SPM_INCDIR}")
  if(ENV_SPM_INCDIR)
    list(APPEND _inc_env "${ENV_SPM_INCDIR}")
  elseif(ENV_SPM_DIR)
    list(APPEND _inc_env "${ENV_SPM_DIR}")
    list(APPEND _inc_env "${ENV_SPM_DIR}/include")
    list(APPEND _inc_env "${ENV_SPM_DIR}/include/spm")
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


  # Try to find the spm header in the given paths
  # ---------------------------------------------------
  # call cmake macro to find the header path
  if(SPM_INCDIR)
    set(SPM_spm.h_DIRS "SPM_spm.h_DIRS-NOTFOUND")
    find_path(SPM_spm.h_DIRS
      NAMES spm.h
      HINTS ${SPM_INCDIR})
  else()
    if(SPM_DIR)
      set(SPM_spm.h_DIRS "SPM_spm.h_DIRS-NOTFOUND")
      find_path(SPM_spm.h_DIRS
        NAMES spm.h
        HINTS ${SPM_DIR}
        PATH_SUFFIXES "include" "include/spm")
    else()
      set(SPM_spm.h_DIRS "SPM_spm.h_DIRS-NOTFOUND")
      find_path(SPM_spm.h_DIRS
        NAMES spm.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "spm")
    endif()
  endif()
  mark_as_advanced(SPM_spm.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (SPM_spm.h_DIRS)
    set(SPM_INCLUDE_DIRS "${SPM_spm.h_DIRS}")
  else ()
    set(SPM_INCLUDE_DIRS "SPM_INCLUDE_DIRS-NOTFOUND")
    if(NOT SPM_FIND_QUIETLY)
      message(STATUS "Looking for spm -- spm.h not found")
    endif()
  endif()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_SPM_LIBDIR "$ENV{SPM_LIBDIR}")
  if(ENV_SPM_LIBDIR)
    list(APPEND _lib_env "${ENV_SPM_LIBDIR}")
  elseif(ENV_SPM_DIR)
    list(APPEND _lib_env "${ENV_SPM_DIR}")
    list(APPEND _lib_env "${ENV_SPM_DIR}/lib")
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

  # Try to find the spm lib in the given paths
  # ------------------------------------------------

  # create list of libs to find
  set(SPM_libs_to_find "spm")
  # call cmake macro to find the lib path
  foreach(spm_lib ${SPM_libs_to_find})
    set(SPM_${spm_lib}_LIBRARY "SPM_${spm_lib}_LIBRARY-NOTFOUND")
    if(SPM_LIBDIR)
      find_library(SPM_${spm_lib}_LIBRARY
        NAMES ${spm_lib}
        HINTS ${SPM_LIBDIR})
    else()
      if(SPM_DIR)
        find_library(SPM_${spm_lib}_LIBRARY
          NAMES ${spm_lib}
          HINTS ${SPM_DIR}
          PATH_SUFFIXES lib lib32 lib64)
      else()
        find_library(SPM_${spm_lib}_LIBRARY
          NAMES ${spm_lib}
          HINTS ${_lib_env})
      endif()
    endif()
  endforeach()
  
  # If found, add path to cmake variable
  # ------------------------------------
  foreach(spm_lib ${SPM_libs_to_find})

    get_filename_component(${spm_lib}_lib_path ${SPM_${spm_lib}_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    if (SPM_LIBRARIES)
      list(APPEND SPM_LIBRARIES "${SPM_${spm_lib}_LIBRARY}")
    else()
      set(SPM_LIBRARIES "${SPM_${spm_lib}_LIBRARY}")
    endif()
    if (SPM_LIBRARY_DIRS)
      list(APPEND SPM_LIBRARY_DIRS "${${spm_lib}_lib_path}")
    else()
      set(SPM_LIBRARY_DIRS "${${spm_lib}_lib_path}")
    endif()
    mark_as_advanced(SPM_${spm_lib}_LIBRARY)

  endforeach(spm_lib ${SPM_libs_to_find})

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT SPM_FOUND) OR (SPM_GIVEN_BY_USER) )


# check a function to validate the find
if(SPM_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # SPM
  if (SPM_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${SPM_INCLUDE_DIRS}")
  endif()
  if (SPM_CFLAGS_OTHER)
    list(APPEND REQUIRED_FLAGS "${SPM_CFLAGS_OTHER}")
  endif()
  if (SPM_LDFLAGS_OTHER)
    list(APPEND REQUIRED_LDFLAGS "${SPM_LDFLAGS_OTHER}")
  endif()
  foreach(libdir ${SPM_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${SPM_LIBRARIES}")
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
  unset(SPM_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(spmInit SPM_WORKS)
  mark_as_advanced(SPM_WORKS)

  if(SPM_WORKS)
    # save link with dependencies
    set(SPM_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(SPM_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(SPM_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(SPM_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(SPM_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT SPM_FIND_QUIETLY)
      message(STATUS "Looking for SPM : test of spm() fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe SPM is linked with specific libraries. "
        "Have you tried with COMPONENTS (MPI/SEQ, PARSEC, STARPU, STARPU_CUDA, SCOTCH, PTSCOTCH, METIS)? "
        "See the explanation in FindSPM.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(SPM_LIBRARIES)


# Check the size of spm_int_t
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${SPM_INCLUDE_DIRS_DEP})

include(CheckCSourceRuns)
set(SPM_C_TEST_SPM_INT_4 "
#include <spm.h>
int main(int argc, char **argv) {
  if (sizeof(spm_int_t) == 4)
    return 0;
  else
    return 1;
}
")
set(SPM_C_TEST_SPM_INT_8 "
#include <spm.h>
int main(int argc, char **argv) {
  if (sizeof(spm_int_t) == 8)
    return 0;
  else
    return 1;
}
")
unset(SPM_INT_4 CACHE)
unset(SPM_INT_8 CACHE)
check_c_source_runs("${SPM_C_TEST_SPM_INT_4}" SPM_INT_4)
check_c_source_runs("${SPM_C_TEST_SPM_INT_8}" SPM_INT_8)
if(NOT SPM_INT_4)
  if(NOT SPM_INT_8)
    set(SPM_INTSIZE -1)
  else()
    set(SPM_INTSIZE 8)
  endif()
else()
  set(SPM_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")

if (SPM_LIBRARIES)
  list(GET SPM_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT SPM_LIBRARY_DIRS)
    set(SPM_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(SPM_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of SPM library" FORCE)
  else()
    set(SPM_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of SPM library" FORCE)
  endif()
endif()
mark_as_advanced(SPM_DIR)
mark_as_advanced(SPM_DIR_FOUND)

if (SPM_LOOK_FOR_FORTRAN)
  find_library(SPM_spmf_LIBRARY
               NAMES spmf
               HINTS ${SPM_LIBRARY_DIRS})
  mark_as_advanced(SPM_spmf_LIBRARY)
  if (SPM_spmf_LIBRARY)
    set(SPM_LIBRARIES_FORTRAN ${SPM_spmf_LIBRARY})
    list(INSERT SPM_LIBRARIES 0 "${SPM_LIBRARIES_FORTRAN}")
    list(INSERT SPM_LIBRARIES_DEP 0 "${SPM_LIBRARIES_FORTRAN}")
  else()
    set(SPM_LIBRARIES_FORTRAN "SPM_LIBRARIES_FORTRAN-NOTFOUND")
  endif()
endif()

# check that SPM has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if (SPM_LOOK_FOR_FORTRAN)
  find_package_handle_standard_args(SPM DEFAULT_MSG
    SPM_LIBRARIES_FORTRAN
    SPM_LIBRARIES
    SPM_WORKS)
else()
  find_package_handle_standard_args(SPM DEFAULT_MSG
    SPM_LIBRARIES
    SPM_WORKS)
endif()
