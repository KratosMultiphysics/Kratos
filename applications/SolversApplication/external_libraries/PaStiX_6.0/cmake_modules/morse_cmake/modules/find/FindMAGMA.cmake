###
# WARNING: not maintained anymore
###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find MAGMA include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(MAGMA
#               [REQUIRED]             # Fail with error if magma is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  MAGMA depends on the following libraries:
#   - CUDA/cuBLAS
#   - LAPACK
#   - CBLAS
#
#  COMPONENTS are optional libraries MAGMA could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - no components are available for now: maybe PLASMA in the future?
#
# Results are reported in variables:
#  MAGMA_FOUND            - True if headers and requested libraries were found
#  MAGMA_C_FLAGS          - list of required compilation flags (excluding -I)
#  MAGMA_LINKER_FLAGS     - list of required linker flags (excluding -l and -L)
#  MAGMA_INCLUDE_DIRS     - magma include directories
#  MAGMA_LIBRARY_DIRS     - Link directories for magma libraries
#  MAGMA_LIBRARIES        - magma libraries
#  MAGMA_INCLUDE_DIRS_DEP - magma + dependencies include directories
#  MAGMA_LIBRARY_DIRS_DEP - magma + dependencies link directories
#  MAGMA_LIBRARIES_DEP    - magma libraries + dependencies
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DMAGMA_DIR=path/to/magma):
#  MAGMA_DIR              - Where to find the base directory of magma
#  MAGMA_INCDIR           - Where to find the header files
#  MAGMA_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: MAGMA_DIR, MAGMA_INCDIR, MAGMA_LIBDIR
#
#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2016 Florent Pruvost
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

if(NOT MAGMA_FOUND)
  set(MAGMA_DIR "" CACHE PATH "Installation directory of MAGMA library")
  if (NOT MAGMA_FIND_QUIETLY)
    message(STATUS "A cache variable, namely MAGMA_DIR, has been set to specify the install directory of MAGMA")
  endif()
endif(NOT MAGMA_FOUND)

# MAGMA depends on CUDA anyway, try to find it
if(MAGMA_FIND_REQUIRED)
  find_package(CUDA REQUIRED)
else()
  find_package(CUDA)
endif()
if( CUDA_FOUND )
  libraries_absolute_path(CUDA_LIBRARIES "$ENV{CUDA_ROOT}/lib64")
endif ()
# MAGMA depends on cuBLAS which should come with CUDA, if not found -> error
if (NOT CUDA_CUBLAS_LIBRARIES)
  if(MAGMA_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for MAGMA - MAGMA depends on cuBLAS which has "
      "not been found (should come with cuda install)")
  endif()
endif()
# MAGMA depends on LAPACK anyway, try to find it
if(MAGMA_FIND_REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()
# MAGMA depends on CBLAS anyway, try to find it
if(MAGMA_FIND_REQUIRED)
  find_package(CBLAS REQUIRED)
else()
  find_package(CBLAS)
endif()

set(ENV_MAGMA_DIR "$ENV{MAGMA_DIR}")
set(ENV_MAGMA_INCDIR "$ENV{MAGMA_INCDIR}")
set(ENV_MAGMA_LIBDIR "$ENV{MAGMA_LIBDIR}")
set(MAGMA_GIVEN_BY_USER "FALSE")
if ( MAGMA_DIR OR ( MAGMA_INCDIR AND MAGMA_LIBDIR) OR ENV_MAGMA_DIR OR (ENV_MAGMA_INCDIR AND ENV_MAGMA_LIBDIR) )
  set(MAGMA_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT MAGMA_GIVEN_BY_USER)

  pkg_search_module(MAGMA magma)
  if (NOT MAGMA_FIND_QUIETLY)
    if (MAGMA_FOUND AND MAGMA_LIBRARIES)
      message(STATUS "Looking for MAGMA - found using PkgConfig")
      #if(NOT MAGMA_INCLUDE_DIRS)
      #    message("${Magenta}MAGMA_INCLUDE_DIRS is empty using PkgConfig."
      #        "Perhaps the path to magma headers is already present in your"
      #        "C(PLUS)_INCLUDE_PATH environment variable.${ColourReset}")
      #endif()
    else()
      message(STATUS "${Magenta}Looking for MAGMA - not found using PkgConfig. "
        "\n   Perhaps you should add the directory containing magma.pc "
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (MAGMA_FIND_VERSION_EXACT)
    if( NOT (MAGMA_FIND_VERSION_MAJOR STREQUAL MAGMA_VERSION_MAJOR) OR
        NOT (MAGMA_FIND_VERSION_MINOR STREQUAL MAGMA_VERSION_MINOR) )
      if(NOT MAGMA_FIND_QUIETLY)
        message(FATAL_ERROR
          "MAGMA version found is ${MAGMA_VERSION_STRING} "
          "when required is ${MAGMA_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (MAGMA_FIND_VERSION_MAJOR STRGREATER MAGMA_VERSION_MAJOR) OR
        (MAGMA_FIND_VERSION_MINOR STRGREATER MAGMA_VERSION_MINOR) )
      if(NOT MAGMA_FIND_QUIETLY)
        message(FATAL_ERROR
          "MAGMA version found is ${MAGMA_VERSION_STRING} "
          "when required is ${MAGMA_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  # if pkg-config is used: these variables are empty
  # the pkg_search_module call will set the following:
  # MAGMA_LDFLAGS: all required linker flags
  # MAGMA_CFLAGS:  all required cflags
  set(MAGMA_INCLUDE_DIRS_DEP "")
  set(MAGMA_LIBRARY_DIRS_DEP "")
  set(MAGMA_LIBRARIES_DEP "")
  # replace it anyway: we should update it with dependencies given by pkg-config
  set(MAGMA_INCLUDE_DIRS_DEP "${MAGMA_STATIC_INCLUDE_DIRS}")
  set(MAGMA_LIBRARY_DIRS_DEP "${MAGMA_STATIC_LIBRARY_DIRS}")
  set(MAGMA_LIBRARIES_DEP "${MAGMA_STATIC_LIBRARIES}")
  set(MAGMA_C_FLAGS "${MAGMA_CFLAGS_OTHER}")

endif(PKG_CONFIG_EXECUTABLE AND NOT MAGMA_GIVEN_BY_USER)

# if MAGMA is not found using pkg-config
if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT MAGMA_FOUND) OR (MAGMA_GIVEN_BY_USER) )

  if (NOT MAGMA_FIND_QUIETLY)
    message(STATUS "Looking for MAGMA - PkgConfig not used")
  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_MAGMA_DIR "$ENV{MAGMA_DIR}")
  set(ENV_MAGMA_INCDIR "$ENV{MAGMA_INCDIR}")
  if(ENV_MAGMA_INCDIR)
    list(APPEND _inc_env "${ENV_MAGMA_INCDIR}")
  elseif(ENV_MAGMA_DIR)
    list(APPEND _inc_env "${ENV_MAGMA_DIR}")
    list(APPEND _inc_env "${ENV_MAGMA_DIR}/include")
    list(APPEND _inc_env "${ENV_MAGMA_DIR}/include/magma")
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


  # Try to find the magma header in the given paths
  # -------------------------------------------------
  # call cmake macro to find the header path
  if(MAGMA_INCDIR)
    set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
    find_path(MAGMA_magma.h_DIRS
      NAMES magma.h
      HINTS ${MAGMA_INCDIR})
  else()
    if(MAGMA_DIR)
      set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
      find_path(MAGMA_magma.h_DIRS
        NAMES magma.h
        HINTS ${MAGMA_DIR}
        PATH_SUFFIXES "include" "include/magma")
    else()
      set(MAGMA_magma.h_DIRS "MAGMA_magma.h_DIRS-NOTFOUND")
      find_path(MAGMA_magma.h_DIRS
        NAMES magma.h
        HINTS ${_inc_env})
    endif()
  endif()
  mark_as_advanced(MAGMA_magma.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (MAGMA_magma.h_DIRS)
    set(MAGMA_INCLUDE_DIRS "${MAGMA_magma.h_DIRS}")
  else ()
    set(MAGMA_INCLUDE_DIRS "MAGMA_INCLUDE_DIRS-NOTFOUND")
    if(NOT MAGMA_FIND_QUIETLY)
      message(STATUS "Looking for magma -- magma.h not found")
    endif()
  endif()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_MAGMA_LIBDIR "$ENV{MAGMA_LIBDIR}")
  if(ENV_MAGMA_LIBDIR)
    list(APPEND _lib_env "${ENV_MAGMA_LIBDIR}")
  elseif(ENV_MAGMA_DIR)
    list(APPEND _lib_env "${ENV_MAGMA_DIR}")
    list(APPEND _lib_env "${ENV_MAGMA_DIR}/lib")
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

  # Try to find the magma lib in the given paths
  # ----------------------------------------------

  # call cmake macro to find the lib path
  if(MAGMA_LIBDIR)
    set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
    find_library(MAGMA_magma_LIBRARY
      NAMES magma
      HINTS ${MAGMA_LIBDIR})
  else()
    if(MAGMA_DIR)
      set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
      find_library(MAGMA_magma_LIBRARY
        NAMES magma
        HINTS ${MAGMA_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    else()
      set(MAGMA_magma_LIBRARY "MAGMA_magma_LIBRARY-NOTFOUND")
      find_library(MAGMA_magma_LIBRARY
        NAMES magma
        HINTS ${_lib_env})
    endif()
  endif()
  mark_as_advanced(MAGMA_magma_LIBRARY)

  # If found, add path to cmake variable
  # ------------------------------------
  if (MAGMA_magma_LIBRARY)
    get_filename_component(magma_lib_path "${MAGMA_magma_LIBRARY}" PATH)
    # set cmake variables
    set(MAGMA_LIBRARIES    "${MAGMA_magma_LIBRARY}")
    set(MAGMA_LIBRARY_DIRS "${magma_lib_path}")
  else ()
    set(MAGMA_LIBRARIES    "MAGMA_LIBRARIES-NOTFOUND")
    set(MAGMA_LIBRARY_DIRS "MAGMA_LIBRARY_DIRS-NOTFOUND")
    if(NOT MAGMA_FIND_QUIETLY)
      message(STATUS "Looking for magma -- lib magma not found")
    endif()
  endif ()

  # check a function to validate the find
  if (MAGMA_LIBRARIES)

    set(REQUIRED_LDFLAGS)
    set(REQUIRED_INCDIRS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # MAGMA
    if (MAGMA_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${MAGMA_INCLUDE_DIRS}")
    endif()
    if (MAGMA_LIBRARY_DIRS)
      set(REQUIRED_LIBDIRS "${MAGMA_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${MAGMA_LIBRARIES}")
    # CBLAS
    if (CBLAS_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS_DEP}")
    elseif (CBLAS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
    endif()
    if(CBLAS_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS_DEP}")
    elseif(CBLAS_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
    endif()
    if (CBLAS_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES_DEP}")
    elseif(CBLAS_LIBRARIES)
      list(APPEND REQUIRED_LIBS "${CBLAS_LIBRARIES}")
    endif()
    if (BLAS_LINKER_FLAGS)
      list(APPEND REQUIRED_LDFLAGS "${BLAS_LINKER_FLAGS}")
    endif()
    # LAPACK
    if (LAPACK_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS}")
    endif()
    if(LAPACK_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES}")
    if (LAPACK_LINKER_FLAGS)
      list(APPEND REQUIRED_LDFLAGS "${LAPACK_LINKER_FLAGS}")
    endif()
    # CUDA
    if (CUDA_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
    endif()
    if(CUDA_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${CUDA_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")
    finds_remove_duplicates()
    # set required libraries for link
    set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
    set(CMAKE_REQUIRED_LIBRARIES)
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LDFLAGS}")
    foreach(lib_dir ${REQUIRED_LIBDIRS})
      list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
    endforeach()
    list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
    string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

    # test link
    unset(MAGMA_WORKS CACHE)
    include(CheckFunctionExists)
    check_function_exists(magma_dgetrf MAGMA_WORKS)
    mark_as_advanced(MAGMA_WORKS)

    if(MAGMA_WORKS)
      # save link with dependencies
      set(MAGMA_LIBRARIES_DEP    "${REQUIRED_LIBS}")
      set(MAGMA_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(MAGMA_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(MAGMA_LINKER_FLAGS     "${REQUIRED_LDFLAGS}")
    else()
      if(NOT MAGMA_FIND_QUIETLY)
        message(STATUS "Looking for magma : test of magma_dgetrf with
                magma, cblas, cuda and lapack libraries fails")
        message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(MAGMA_LIBRARIES)

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT MAGMA_FOUND) OR (MAGMA_GIVEN_BY_USER) )

if (MAGMA_LIBRARIES)
  if (MAGMA_LIBRARY_DIRS)
    foreach(dir ${MAGMA_LIBRARY_DIRS})
      if ("${dir}" MATCHES "magma")
        set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET MAGMA_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(MAGMA_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(MAGMA_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of MAGMA library" FORCE)
  else()
    set(MAGMA_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of MAGMA library" FORCE)
  endif()
endif()
mark_as_advanced(MAGMA_DIR)
mark_as_advanced(MAGMA_DIR_FOUND)

# check that MAGMA has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
if (PKG_CONFIG_EXECUTABLE AND MAGMA_FOUND)
  find_package_handle_standard_args(MAGMA DEFAULT_MSG
    MAGMA_LIBRARIES)
else()
  find_package_handle_standard_args(MAGMA DEFAULT_MSG
    MAGMA_LIBRARIES
    MAGMA_WORKS)
endif()
