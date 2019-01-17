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
# - Find CBLAS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(CBLAS
#               [REQUIRED] # Fail with error if cblas is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  CBLAS depends on the following libraries:
#   - BLAS
#
# This module finds headers and cblas library.
# Results are reported in variables:
#  CBLAS_FOUND             - True if headers and requested libraries were found
#  CBLAS_CFLAGS_OTHER      - cblas compiler flags without headers paths
#  CBLAS_LDFLAGS_OTHER     - cblas linker flags without libraries
#  CBLAS_INCLUDE_DIRS      - cblas include directories
#  CBLAS_LIBRARY_DIRS      - cblas link directories
#  CBLAS_LIBRARIES         - cblas libraries to be linked (absolute path)
#  CBLAS_CFLAGS_OTHER_DEP  - cblas + dependencies compiler flags without headers paths
#  CBLAS_LDFLAGS_OTHER_DEP - cblas + dependencies linker flags without libraries
#  CBLAS_INCLUDE_DIRS_DEP  - cblas + dependencies include directories
#  CBLAS_LIBRARY_DIRS_DEP  - cblas + dependencies link directories
#  CBLAS_LIBRARIES_DEP     - cblas + dependencies libraries
#  CBLAS_HAS_ZGEMM3M       - True if cblas contains zgemm3m fast complex mat-mat product
#
#  CBLAS_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = CBLAS
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
# options at configure (ex: cmake path/to/project -DCBLAS_DIR=path/to/cblas):
#  CBLAS_DIR              - Where to find the base directory of cblas
#  CBLAS_INCDIR           - Where to find the header files
#  CBLAS_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: CBLAS_DIR, CBLAS_INCDIR, CBLAS_LIBDIR
#
# CBLAS could be directly embedded in BLAS library (ex: Intel MKL) so that
# we test a cblas function with the blas libraries found and set CBLAS
# variables to BLAS ones if test is successful. To skip this feature and
# look for a stand alone cblas, please add the following in your
# CMakeLists.txt before to call find_package(CBLAS):
# set(CBLAS_STANDALONE TRUE)
###
# We handle different modes to find the dependency
#
# - Detection if already installed on the system
#   - CBLAS libraries can be detected from different ways
#     Here is the order of precedence:
#     1) we look in cmake variable CBLAS_LIBDIR or CBLAS_DIR (we guess the libdirs) if defined
#     2) we look in environment variable CBLAS_LIBDIR or CBLAS_DIR (we guess the libdirs) if defined
#     3) we look in common environnment variables depending on the system (INCLUDE, C_INCLUDE_PATH, CPATH - LIB, DYLD_LIBRARY_PATH, LD_LIBRARY_PATH)
#     4) we look in common system paths depending on the system, see for example paths contained in the following cmake variables:
#       - CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES, CMAKE_C_IMPLICIT_LINK_DIRECTORIES
#

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

if (NOT CBLAS_FOUND)
  set(CBLAS_DIR "" CACHE PATH "Installation directory of CBLAS library")
  if (NOT CBLAS_FIND_QUIETLY)
    message(STATUS "A cache variable, namely CBLAS_DIR, has been set to specify the install directory of CBLAS")
  endif()
endif()


# CBLAS depends on BLAS anyway, try to find it
if(CBLAS_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
else()
  find_package(BLAS)
endif()


# find CBLAS
if (BLAS_FOUND)

  if (NOT CBLAS_STANDALONE)
    # check if a cblas function exists in the BLAS lib
    # this can be the case with libs such as MKL, ACML
    include(CheckFunctionExists)
    set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LDFLAGS_OTHER};${BLAS_LIBRARIES}")
    set(CMAKE_REQUIRED_FLAGS "${BLAS_CFLAGS_OTHER}")
    set(CMAKE_REQUIRED_INCLUDES "${BLAS_INCLUDE_DIRS}")
    unset(CBLAS_WORKS CACHE)
    check_function_exists(cblas_dscal CBLAS_WORKS)
    check_function_exists(cblas_zgemm3m CBLAS_ZGEMM3M_FOUND)
    mark_as_advanced(CBLAS_WORKS)
    set(CMAKE_REQUIRED_LIBRARIES)

    if(CBLAS_WORKS)

      # Check for faster complex GEMM routine
      # (only C/Z, no S/D version)
      if ( CBLAS_ZGEMM3M_FOUND )
        add_definitions(-DCBLAS_HAS_ZGEMM3M -DCBLAS_HAS_CGEMM3M)
      endif()

      if(NOT CBLAS_FIND_QUIETLY)
        message(STATUS "Looking for cblas: test with blas succeeds")
      endif()
      # test succeeds: CBLAS is in BLAS
      set(CBLAS_LIBRARIES "${BLAS_LIBRARIES}")

    endif()

  endif (NOT CBLAS_STANDALONE)

  # test fails with blas: try to find CBLAS lib exterior to BLAS
  if (CBLAS_STANDALONE OR NOT CBLAS_WORKS)

    if(NOT CBLAS_WORKS AND NOT CBLAS_FIND_QUIETLY)
      message(STATUS "Looking for cblas : test with blas fails")
    endif()

    # try with pkg-config
    set(ENV_CBLAS_DIR "$ENV{CBLAS_DIR}")
    set(ENV_MKL_DIR "$ENV{MKLROOT}")
    set(ENV_CBLAS_INCDIR "$ENV{CBLAS_INCDIR}")
    set(ENV_CBLAS_LIBDIR "$ENV{CBLAS_LIBDIR}")
    set(CBLAS_GIVEN_BY_USER "FALSE")
    if ( CBLAS_DIR OR ( CBLAS_INCDIR AND CBLAS_LIBDIR) OR ENV_CBLAS_DIR OR ENV_MKL_DIR OR (ENV_CBLAS_INCDIR AND ENV_CBLAS_LIBDIR) )
      set(CBLAS_GIVEN_BY_USER "TRUE")
    endif()

    include(FindPkgConfig)
    find_package(PkgConfig QUIET)
    if( PKG_CONFIG_EXECUTABLE AND NOT CBLAS_GIVEN_BY_USER)

      if (BLA_STATIC)
        set(MKL_STR_BLA_STATIC "static")
      else()
        set(MKL_STR_BLA_STATIC "dynamic")
      endif()
      # try different blas
      if (BLA_VENDOR STREQUAL "Intel10_64lp")
        pkg_search_module(CBLAS mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
      elseif(BLA_VENDOR STREQUAL "Intel10_64lp_seq")
        pkg_search_module(CBLAS mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      elseif(BLA_VENDOR STREQUAL "Open")
        pkg_search_module(CBLAS openblas)
      elseif(BLA_VENDOR STREQUAL "Generic")
        pkg_search_module(CBLAS cblas)
      else()
        pkg_search_module(CBLAS cblas)
        pkg_search_module(CBLAS openblas)
        pkg_search_module(CBLAS mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      endif()

      if (NOT CBLAS_FIND_QUIETLY)
        if (CBLAS_FOUND AND CBLAS_LIBRARIES)
          message(STATUS "Looking for CBLAS - found using PkgConfig")
        else()
          message(STATUS "${Magenta}Looking for CBLAS - not found using PkgConfig."
            "\n   Perhaps you should add the directory containing cblas.pc to"
            "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
        endif()
      endif()

      if (CBLAS_FOUND AND CBLAS_LIBRARIES)
        set(CBLAS_FOUND_WITH_PKGCONFIG "TRUE")
        find_pkgconfig_libraries_absolute_path(CBLAS)
      else()
        set(CBLAS_FOUND_WITH_PKGCONFIG "FALSE")
      endif()

    endif()

    if (NOT CBLAS_FOUND_WITH_PKGCONFIG OR CBLAS_GIVEN_BY_USER)
      # Try to find CBLAS lib
      #######################

      # Looking for include
      # -------------------

      # Add system include paths to search include
      # ------------------------------------------
      unset(_inc_env)
      set(ENV_CBLAS_DIR "$ENV{CBLAS_DIR}")
      set(ENV_CBLAS_INCDIR "$ENV{CBLAS_INCDIR}")
      if(ENV_CBLAS_INCDIR)
        list(APPEND _inc_env "${ENV_CBLAS_INCDIR}")
      elseif(ENV_CBLAS_DIR)
        list(APPEND _inc_env "${ENV_CBLAS_DIR}")
        list(APPEND _inc_env "${ENV_CBLAS_DIR}/include")
        list(APPEND _inc_env "${ENV_CBLAS_DIR}/include/cblas")
      else()
        if(WIN32)
          string(REPLACE ":" ";" _path_env "$ENV{INCLUDE}")
          list(APPEND _inc_env "${_path_env}")
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


      # Try to find the cblas header in the given paths
      # -------------------------------------------------
      # call cmake macro to find the header path
      if(CBLAS_INCDIR)
        set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
        find_path(CBLAS_cblas.h_DIRS
          NAMES cblas.h
          HINTS ${CBLAS_INCDIR})
      else()
        if(CBLAS_DIR)
          set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
          find_path(CBLAS_cblas.h_DIRS
            NAMES cblas.h
            HINTS ${CBLAS_DIR}
            PATH_SUFFIXES "include" "include/cblas")
        else()
          set(CBLAS_cblas.h_DIRS "CBLAS_cblas.h_DIRS-NOTFOUND")
          find_path(CBLAS_cblas.h_DIRS
            NAMES cblas.h
            HINTS ${_inc_env}
            PATH_SUFFIXES "cblas")
        endif()
      endif()
      mark_as_advanced(CBLAS_cblas.h_DIRS)

      # If found, add path to cmake variable
      # ------------------------------------
      if (CBLAS_cblas.h_DIRS)
        set(CBLAS_INCLUDE_DIRS "${CBLAS_cblas.h_DIRS}")
      else ()
        set(CBLAS_INCLUDE_DIRS "CBLAS_INCLUDE_DIRS-NOTFOUND")
        if(NOT CBLAS_FIND_QUIETLY)
          message(STATUS "Looking for cblas -- cblas.h not found")
        endif()
      endif()


      # Looking for lib
      # ---------------

      # Add system library paths to search lib
      # --------------------------------------
      unset(_lib_env)
      set(ENV_CBLAS_LIBDIR "$ENV{CBLAS_LIBDIR}")
      if(ENV_CBLAS_LIBDIR)
        list(APPEND _lib_env "${ENV_CBLAS_LIBDIR}")
      elseif(ENV_CBLAS_DIR)
        list(APPEND _lib_env "${ENV_CBLAS_DIR}")
        list(APPEND _lib_env "${ENV_CBLAS_DIR}/lib")
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

      # Try to find the cblas lib in the given paths
      # ----------------------------------------------

      # call cmake macro to find the lib path
      if(CBLAS_LIBDIR)
        set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
        find_library(CBLAS_cblas_LIBRARY
          NAMES cblas
          HINTS ${CBLAS_LIBDIR})
      else()
        if(CBLAS_DIR)
          set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
          find_library(CBLAS_cblas_LIBRARY
            NAMES cblas
            HINTS ${CBLAS_DIR}
            PATH_SUFFIXES lib lib32 lib64)
        else()
          set(CBLAS_cblas_LIBRARY "CBLAS_cblas_LIBRARY-NOTFOUND")
          find_library(CBLAS_cblas_LIBRARY
            NAMES cblas
            HINTS ${_lib_env})
        endif()
      endif()
      mark_as_advanced(CBLAS_cblas_LIBRARY)

      # If found, add path to cmake variable
      # ------------------------------------
      if (CBLAS_cblas_LIBRARY)
        get_filename_component(cblas_lib_path "${CBLAS_cblas_LIBRARY}" PATH)
        # set cmake variables
        set(CBLAS_LIBRARIES    "${CBLAS_cblas_LIBRARY}")
        set(CBLAS_LIBRARY_DIRS "${cblas_lib_path}")
      else ()
        set(CBLAS_LIBRARIES    "CBLAS_LIBRARIES-NOTFOUND")
        set(CBLAS_LIBRARY_DIRS "CBLAS_LIBRARY_DIRS-NOTFOUND")
        if (NOT CBLAS_FIND_QUIETLY)
          message(STATUS "Looking for cblas -- lib cblas not found")
        endif()
      endif ()
    endif (NOT CBLAS_FOUND_WITH_PKGCONFIG OR CBLAS_GIVEN_BY_USER)

  endif (CBLAS_STANDALONE OR NOT CBLAS_WORKS)

  # check a function to validate the find
  if(CBLAS_LIBRARIES)

    set(REQUIRED_INCDIRS)
    set(REQUIRED_FLAGS)
    set(REQUIRED_LDFLAGS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # CBLAS
    if (CBLAS_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${CBLAS_INCLUDE_DIRS}")
    endif()
    if (CBLAS_CFLAGS_OTHER)
      set(REQUIRED_FLAGS "${CBLAS_CFLAGS_OTHER}")
    endif()
    if (CBLAS_LDFLAGS_OTHER)
      set(REQUIRED_LDFLAGS "${CBLAS_LDFLAGS_OTHER}")
    endif()
    if (CBLAS_LIBRARY_DIRS)
      set(REQUIRED_LIBDIRS "${CBLAS_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${CBLAS_LIBRARIES}")
    # BLAS
    if (BLAS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${BLAS_INCLUDE_DIRS}")
    endif()
    if (BLAS_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${BLAS_CFLAGS_OTHER}")
    endif()
    if (BLAS_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${BLAS_LDFLAGS_OTHER}")
    endif()
    if (BLAS_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${BLAS_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${BLAS_LIBRARIES}")

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
    unset(CBLAS_WORKS CACHE)
    include(CheckFunctionExists)
    check_function_exists(cblas_dscal CBLAS_WORKS)
    mark_as_advanced(CBLAS_WORKS)

    if(CBLAS_WORKS)

      # Check for faster complex GEMM routine
      # (only C/Z, no S/D version)
      check_function_exists(cblas_zgemm3m CBLAS_ZGEMM3M_FOUND)
      if ( CBLAS_ZGEMM3M_FOUND )
        add_definitions(-DCBLAS_HAS_ZGEMM3M -DCBLAS_HAS_CGEMM3M)
      endif()

      # save link with dependencies
      set(CBLAS_LIBRARIES_DEP "${REQUIRED_LIBS}")
      set(CBLAS_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(CBLAS_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(CBLAS_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
      set(CBLAS_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
    else()
      if(NOT CBLAS_FIND_QUIETLY)
        message(STATUS "Looking for cblas : test of cblas_dscal with cblas and blas libraries fails")
        message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(CBLAS_LIBRARIES)

else(BLAS_FOUND)

  if (NOT CBLAS_FIND_QUIETLY)
    message(STATUS "CBLAS requires BLAS but BLAS has not been found."
      "Please look for BLAS first.")
  endif()

endif(BLAS_FOUND)

if (CBLAS_LIBRARIES)
  list(GET CBLAS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT CBLAS_LIBRARY_DIRS)
    set(CBLAS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(CBLAS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of CBLAS library" FORCE)
  else()
    set(CBLAS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of CBLAS library" FORCE)
  endif()
endif()
mark_as_advanced(CBLAS_DIR)
mark_as_advanced(CBLAS_DIR_FOUND)

# check that CBLAS has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS DEFAULT_MSG
  CBLAS_LIBRARIES
  CBLAS_WORKS)
