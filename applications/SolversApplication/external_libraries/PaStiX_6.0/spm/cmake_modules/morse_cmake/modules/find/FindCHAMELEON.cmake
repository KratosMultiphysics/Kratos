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
# - Find CHAMELEON include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(CHAMELEON
#               [REQUIRED]             # Fail with error if chameleon is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  CHAMELEON depends on the following libraries:
#   - Threads, m, rt
#   - HWLOC
#   - CBLAS
#   - LAPACKE
#   - TMG
#   - At least one runtime, default is StarPU
#     (For QUARK, use COMPONENTS QUARK)
#
#  COMPONENTS are optional libraries CHAMELEON could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - STARPU (default): to activate detection of Chameleon linked with StarPU
#   - QUARK (STARPU will be deactivated): to activate detection of Chameleon linked with QUARK
#   - CUDA (comes with cuBLAS): to activate detection of Chameleon linked with CUDA
#   - MPI: to activate detection of Chameleon linked with MPI
#   - FXT: to activate detection of Chameleon linked with StarPU+FXT
#
# This module finds headers and chameleon library.
# Results are reported in variables:
#  CHAMELEON_FOUND             - True if headers and requested libraries were found
#  CHAMELEON_CFLAGS_OTHER      - chameleon compiler flags without headers paths
#  CHAMELEON_LDFLAGS_OTHER     - chameleon linker flags without libraries
#  CHAMELEON_INCLUDE_DIRS      - chameleon include directories
#  CHAMELEON_LIBRARY_DIRS      - chameleon link directories
#  CHAMELEON_LIBRARIES         - chameleon libraries to be linked (absolute path)
#  CHAMELEON_CFLAGS_OTHER_DEP  - chameleon + dependencies compiler flags without headers paths
#  CHAMELEON_LDFLAGS_OTHER_DEP - chameleon + dependencies linker flags without libraries
#  CHAMELEON_INCLUDE_DIRS_DEP  - chameleon + dependencies include directories
#  CHAMELEON_LIBRARY_DIRS_DEP  - chameleon + dependencies link directories
#  CHAMELEON_LIBRARIES_DEP     - chameleon + dependencies libraries
#
#  CHAMELEON_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = CHAMELEON
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
# options at configure (ex: cmake path/to/project -DCHAMELEON_DIR=path/to/chameleon):
#  CHAMELEON_DIR              - Where to find the base directory of chameleon
#  CHAMELEON_INCDIR           - Where to find the header files
#  CHAMELEON_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: CHAMELEON_DIR, CHAMELEON_INCDIR, CHAMELEON_LIBDIR

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

if (NOT CHAMELEON_FOUND)
  set(CHAMELEON_DIR "" CACHE PATH "Installation directory of CHAMELEON library")
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "A cache variable, namely CHAMELEON_DIR, has been set to specify the install directory of CHAMELEON")
  endif()
endif()

# Try to find CHAMELEON dependencies if specified as COMPONENTS during the call
set(CHAMELEON_LOOK_FOR_STARPU ON)
set(CHAMELEON_LOOK_FOR_QUARK OFF)
set(CHAMELEON_LOOK_FOR_CUDA OFF)
set(CHAMELEON_LOOK_FOR_MPI OFF)
set(CHAMELEON_LOOK_FOR_FXT OFF)

if( CHAMELEON_FIND_COMPONENTS )
  foreach( component ${CHAMELEON_FIND_COMPONENTS} )
    if (${component} STREQUAL "STARPU")
      # means we look for Chameleon with StarPU
      set(CHAMELEON_LOOK_FOR_STARPU ON)
      set(CHAMELEON_LOOK_FOR_QUARK OFF)
    endif()
    if (${component} STREQUAL "QUARK")
      # means we look for Chameleon with QUARK
      set(CHAMELEON_LOOK_FOR_QUARK ON)
      set(CHAMELEON_LOOK_FOR_STARPU OFF)
    endif()
    if (${component} STREQUAL "CUDA")
      # means we look for Chameleon with CUDA
      set(CHAMELEON_LOOK_FOR_CUDA ON)
    endif()
    if (${component} STREQUAL "MPI")
      # means we look for Chameleon with MPI
      set(CHAMELEON_LOOK_FOR_MPI ON)
    endif()
    if (${component} STREQUAL "FXT")
      # means we look for Chameleon with FXT
      set(CHAMELEON_LOOK_FOR_FXT ON)
    endif()
  endforeach()
endif()

set(ENV_CHAMELEON_DIR "$ENV{CHAMELEON_DIR}")
set(ENV_CHAMELEON_INCDIR "$ENV{CHAMELEON_INCDIR}")
set(ENV_CHAMELEON_LIBDIR "$ENV{CHAMELEON_LIBDIR}")
set(CHAMELEON_GIVEN_BY_USER "FALSE")
if ( CHAMELEON_DIR OR ( CHAMELEON_INCDIR AND CHAMELEON_LIBDIR) OR ENV_CHAMELEON_DIR OR (ENV_CHAMELEON_INCDIR AND ENV_CHAMELEON_LIBDIR) )
  set(CHAMELEON_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT CHAMELEON_GIVEN_BY_USER)

  pkg_search_module(CHAMELEON chameleon)

  if (NOT CHAMELEON_FIND_QUIETLY)
    if (CHAMELEON_FOUND AND CHAMELEON_LIBRARIES)
      message(STATUS "Looking for CHAMELEON - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for CHAMELEON - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing chameleon.pc"
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (CHAMELEON_FIND_VERSION_EXACT)
    if( NOT (CHAMELEON_FIND_VERSION_MAJOR STREQUAL CHAMELEON_VERSION_MAJOR) OR
        NOT (CHAMELEON_FIND_VERSION_MINOR STREQUAL CHAMELEON_VERSION_MINOR) )
      if(NOT CHAMELEON_FIND_QUIETLY)
        message(FATAL_ERROR
          "CHAMELEON version found is ${CHAMELEON_VERSION_STRING}"
          "when required is ${CHAMELEON_FIND_VERSION}")
      endif()
    endif()
  else()
    # if the version found is older than the required then error
    if( (CHAMELEON_FIND_VERSION_MAJOR STRGREATER CHAMELEON_VERSION_MAJOR) OR
        (CHAMELEON_FIND_VERSION_MINOR STRGREATER CHAMELEON_VERSION_MINOR) )
      if(NOT CHAMELEON_FIND_QUIETLY)
        message(FATAL_ERROR
          "CHAMELEON version found is ${CHAMELEON_VERSION_STRING}"
          "when required is ${CHAMELEON_FIND_VERSION} or newer")
      endif()
    endif()
  endif()

  if (CHAMELEON_FOUND AND CHAMELEON_LIBRARIES)
    set(CHAMELEON_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(CHAMELEON)
  else()
    set(CHAMELEON_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT CHAMELEON_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT CHAMELEON_FOUND) OR (CHAMELEON_GIVEN_BY_USER) )

  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - PkgConfig not used")
  endif()

  # Dependencies detection
  # ----------------------

  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect pthread")
  endif()
  if (CHAMELEON_FIND_REQUIRED)
    find_package(Threads REQUIRED)
  else()
    find_package(Threads)
  endif()
  set(CHAMELEON_EXTRA_LIBRARIES "")
  if( THREADS_FOUND )
    libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
    list(APPEND CHAMELEON_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
  endif ()

  # Add math library to the list of extra
  # it normally exists on all common systems provided with a C compiler
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect libm")
  endif()
  set(CHAMELEON_M_LIBRARIES "")
  if(UNIX OR WIN32)
    find_library(
      CHAMELEON_M_m_LIBRARY
      NAMES m
      )
    mark_as_advanced(CHAMELEON_M_m_LIBRARY)
    if (CHAMELEON_M_m_LIBRARY)
      list(APPEND CHAMELEON_M_LIBRARIES "${CHAMELEON_M_m_LIBRARY}")
      list(APPEND CHAMELEON_EXTRA_LIBRARIES "${CHAMELEON_M_m_LIBRARY}")
    else()
      if (CHAMELEON_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find libm on your system."
          "Are you sure to a have a C compiler installed?")
      endif()
    endif()
  endif()

  # Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
  # on Unix systems except Apple ones because it does not exist on it
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect librt")
  endif()
  set(CHAMELEON_RT_LIBRARIES "")
  if(UNIX AND NOT APPLE)
    find_library(
      CHAMELEON_RT_rt_LIBRARY
      NAMES rt
      )
    mark_as_advanced(CHAMELEON_RT_rt_LIBRARY)
    if (CHAMELEON_RT_rt_LIBRARY)
      list(APPEND CHAMELEON_RT_LIBRARIES "${CHAMELEON_RT_rt_LIBRARY}")
      list(APPEND CHAMELEON_EXTRA_LIBRARIES "${CHAMELEON_RT_rt_LIBRARY}")
    else()
      if (CHAMELEON_FIND_REQUIRED)
        message(FATAL_ERROR "Could NOT find librt on your system")
      endif()
    endif()
  endif()

  # CHAMELEON depends on CBLAS
  #---------------------------
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect CBLAS (depends on BLAS)")
  endif()
  if (CHAMELEON_FIND_REQUIRED)
    find_package(CBLAS REQUIRED)
  else()
    find_package(CBLAS)
  endif()

  # CHAMELEON depends on LAPACKE
  #-----------------------------

  # standalone version of lapacke seems useless for now
  # let the comment in case we meet some problems of non existing lapacke
  # functions in lapack library such as mkl, acml, ...
  #set(LAPACKE_STANDALONE TRUE)
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect LAPACKE (depends on LAPACK)")
  endif()
  if (CHAMELEON_FIND_REQUIRED)
    find_package(LAPACKE REQUIRED)
  else()
    find_package(LAPACKE)
  endif()

  # CHAMELEON depends on TMG
  #-------------------------
  if (NOT CHAMELEON_FIND_QUIETLY)
    message(STATUS "Looking for CHAMELEON - Try to detect TMG (depends on LAPACK)")
  endif()
  if (CHAMELEON_FIND_REQUIRED)
    find_package(TMG REQUIRED)
  else()
    find_package(TMG)
  endif()

  # CHAMELEON may depend on CUDA/CUBLAS
  #------------------------------------
  if (NOT CUDA_FOUND AND CHAMELEON_LOOK_FOR_CUDA)
    if (CHAMELEON_FIND_REQUIRED AND CHAMELEON_FIND_REQUIRED_CUDA)
      find_package(CUDA REQUIRED)
    else()
      find_package(CUDA)
    endif()
    if (CUDA_FOUND)
      mark_as_advanced(CUDA_BUILD_CUBIN)
      mark_as_advanced(CUDA_BUILD_EMULATION)
      mark_as_advanced(CUDA_SDK_ROOT_DIR)
      mark_as_advanced(CUDA_TOOLKIT_ROOT_DIR)
      mark_as_advanced(CUDA_VERBOSE_BUILD)
      libraries_absolute_path(CUDA_LIBRARIES "$ENV{CUDA_ROOT}/lib64")
    endif()
  endif()

  # CHAMELEON depends on MPI
  #-------------------------
  if( NOT MPI_FOUND AND CHAMELEON_LOOK_FOR_MPI )

    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if(NOT MPI_C_COMPILER)
      set(MPI_C_COMPILER mpicc)
    endif()
    if (CHAMELEON_FIND_REQUIRED AND CHAMELEON_FIND_REQUIRED_MPI)
      find_package(MPI REQUIRED)
    else()
      find_package(MPI)
    endif()
    if (MPI_FOUND)
      mark_as_advanced(MPI_LIBRARY)
      mark_as_advanced(MPI_EXTRA_LIBRARY)
    endif()

  endif()

  if( NOT STARPU_FOUND AND CHAMELEON_LOOK_FOR_STARPU )

    set(CHAMELEON_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

    # create list of components in order to make a single call to find_package(starpu...)
    # we explicitly need a StarPU version built with hwloc
    set(STARPU_COMPONENT_LIST "HWLOC")

    # StarPU may depend on MPI
    # allows to use an external mpi compilation by setting compilers with
    # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
    # at cmake configure
    if (CHAMELEON_LOOK_FOR_MPI)
      if(NOT MPI_C_COMPILER)
        set(MPI_C_COMPILER mpicc)
      endif()
      list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if (CHAMELEON_LOOK_FOR_CUDA)
      list(APPEND STARPU_COMPONENT_LIST "CUDA")
    endif()
    if (CHAMELEON_LOOK_FOR_FXT)
      list(APPEND STARPU_COMPONENT_LIST "FXT")
    endif()
    if (CHAMELEON_FIND_REQUIRED AND CHAMELEON_FIND_REQUIRED_STARPU)
      find_package(STARPU ${CHAMELEON_STARPU_VERSION} REQUIRED
        COMPONENTS ${STARPU_COMPONENT_LIST})
    else()
      find_package(STARPU ${CHAMELEON_STARPU_VERSION}
        COMPONENTS ${STARPU_COMPONENT_LIST})
    endif()

  endif()

  if( NOT QUARK_FOUND AND CHAMELEON_LOOK_FOR_QUARK )

    # try to find quark runtime
    if (CHAMELEON_FIND_REQUIRED AND CHAMELEON_FIND_REQUIRED_QUARK)
      find_package(QUARK REQUIRED COMPONENTS HWLOC)
    else()
      find_package(QUARK COMPONENTS HWLOC)
    endif()

  endif()

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_CHAMELEON_DIR "$ENV{CHAMELEON_DIR}")
  set(ENV_CHAMELEON_INCDIR "$ENV{CHAMELEON_INCDIR}")
  if(ENV_CHAMELEON_INCDIR)
    list(APPEND _inc_env "${ENV_CHAMELEON_INCDIR}")
  elseif(ENV_CHAMELEON_DIR)
    list(APPEND _inc_env "${ENV_CHAMELEON_DIR}")
    list(APPEND _inc_env "${ENV_CHAMELEON_DIR}/include")
    list(APPEND _inc_env "${ENV_CHAMELEON_DIR}/include/chameleon")
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


  # Try to find the chameleon header in the given paths
  # ---------------------------------------------------
  # call cmake macro to find the header path
  if(CHAMELEON_INCDIR)
    set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
    find_path(CHAMELEON_morse.h_DIRS
      NAMES morse.h
      HINTS ${CHAMELEON_INCDIR})
  else()
    if(CHAMELEON_DIR)
      set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
      find_path(CHAMELEON_morse.h_DIRS
        NAMES morse.h
        HINTS ${CHAMELEON_DIR}
        PATH_SUFFIXES "include" "include/chameleon")
    else()
      set(CHAMELEON_morse.h_DIRS "CHAMELEON_morse.h_DIRS-NOTFOUND")
      find_path(CHAMELEON_morse.h_DIRS
        NAMES morse.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "chameleon")
    endif()
  endif()
  mark_as_advanced(CHAMELEON_morse.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (CHAMELEON_morse.h_DIRS)
    set(CHAMELEON_INCLUDE_DIRS "${CHAMELEON_morse.h_DIRS}")
  else ()
    set(CHAMELEON_INCLUDE_DIRS "CHAMELEON_INCLUDE_DIRS-NOTFOUND")
    if(NOT CHAMELEON_FIND_QUIETLY)
      message(STATUS "Looking for chameleon -- morse.h not found")
    endif()
  endif()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_CHAMELEON_LIBDIR "$ENV{CHAMELEON_LIBDIR}")
  if(ENV_CHAMELEON_LIBDIR)
    list(APPEND _lib_env "${ENV_CHAMELEON_LIBDIR}")
  elseif(ENV_CHAMELEON_DIR)
    list(APPEND _lib_env "${ENV_CHAMELEON_DIR}")
    list(APPEND _lib_env "${ENV_CHAMELEON_DIR}/lib")
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

  # Try to find the chameleon lib in the given paths
  # ------------------------------------------------

  # create list of libs to find
  set(CHAMELEON_libs_to_find "chameleon")
  if (STARPU_FOUND)
    list(APPEND CHAMELEON_libs_to_find "chameleon_starpu")
  elseif (QUARK_FOUND)
    list(APPEND CHAMELEON_libs_to_find "chameleon_quark")
  endif()
  list(APPEND CHAMELEON_libs_to_find "coreblas")
  if (CHAMELEON_LOOK_FOR_CUDA)
    list(APPEND CHAMELEON_libs_to_find "cudablas")
  endif()

  # call cmake macro to find the lib path
  if(CHAMELEON_LIBDIR)
    foreach(chameleon_lib ${CHAMELEON_libs_to_find})
      set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
      find_library(CHAMELEON_${chameleon_lib}_LIBRARY
        NAMES ${chameleon_lib}
        HINTS ${CHAMELEON_LIBDIR})
    endforeach()
  else()
    if(CHAMELEON_DIR)
      foreach(chameleon_lib ${CHAMELEON_libs_to_find})
        set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
        find_library(CHAMELEON_${chameleon_lib}_LIBRARY
          NAMES ${chameleon_lib}
          HINTS ${CHAMELEON_DIR}
          PATH_SUFFIXES lib lib32 lib64)
      endforeach()
    else()
      foreach(chameleon_lib ${CHAMELEON_libs_to_find})
        set(CHAMELEON_${chameleon_lib}_LIBRARY "CHAMELEON_${chameleon_lib}_LIBRARY-NOTFOUND")
        find_library(CHAMELEON_${chameleon_lib}_LIBRARY
          NAMES ${chameleon_lib}
          HINTS ${_lib_env})
      endforeach()
    endif()
  endif()

  # If found, add path to cmake variable
  # ------------------------------------
  foreach(chameleon_lib ${CHAMELEON_libs_to_find})

    get_filename_component(${chameleon_lib}_lib_path ${CHAMELEON_${chameleon_lib}_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    if (CHAMELEON_LIBRARIES)
      list(APPEND CHAMELEON_LIBRARIES "${CHAMELEON_${chameleon_lib}_LIBRARY}")
    else()
      set(CHAMELEON_LIBRARIES "${CHAMELEON_${chameleon_lib}_LIBRARY}")
    endif()
    if (CHAMELEON_LIBRARY_DIRS)
      list(APPEND CHAMELEON_LIBRARY_DIRS "${${chameleon_lib}_lib_path}")
    else()
      set(CHAMELEON_LIBRARY_DIRS "${${chameleon_lib}_lib_path}")
    endif()
    mark_as_advanced(CHAMELEON_${chameleon_lib}_LIBRARY)

  endforeach(chameleon_lib ${CHAMELEON_libs_to_find})

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT CHAMELEON_FOUND) OR (CHAMELEON_GIVEN_BY_USER) )

# check a function to validate the find
if(CHAMELEON_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # CHAMELEON
  if (CHAMELEON_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${CHAMELEON_INCLUDE_DIRS}")
  endif()
  if (CHAMELEON_CFLAGS_OTHER)
    list(APPEND REQUIRED_FLAGS "${CHAMELEON_CFLAGS_OTHER}")
  endif()
  if (CHAMELEON_LDFLAGS_OTHER)
    list(APPEND REQUIRED_LDFLAGS "${CHAMELEON_LDFLAGS_OTHER}")
  endif()
  foreach(libdir ${CHAMELEON_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${CHAMELEON_LIBRARIES}")
  # STARPU
  if (STARPU_FOUND AND CHAMELEON_LOOK_FOR_STARPU)
    if (STARPU_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${STARPU_INCLUDE_DIRS_DEP}")
    endif()
    if (STARPU_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${STARPU_CFLAGS_OTHER_DEP}")
    endif()
    if (STARPU_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${STARPU_LDFLAGS_OTHER_DEP}")
    endif()
    if(STARPU_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${STARPU_LIBRARY_DIRS_DEP}")
    endif()
    if (STARPU_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${STARPU_LIBRARIES_DEP}")
    endif()
  endif()
  # QUARK
  if (QUARK_FOUND AND CHAMELEON_LOOK_FOR_QUARK)
    if (QUARK_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${QUARK_INCLUDE_DIRS_DEP}")
    endif()
    if (QUARK_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${QUARK_CFLAGS_OTHER_DEP}")
    endif()
    if (QUARK_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${QUARK_LDFLAGS_OTHER_DEP}")
    endif()
    if(QUARK_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${QUARK_LIBRARY_DIRS_DEP}")
    endif()
    if (QUARK_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBS "${QUARK_LIBRARIES_DEP}")
    endif()
  endif()
  # CUDA
  if (CUDA_FOUND AND CHAMELEON_LOOK_FOR_CUDA)
    if (CUDA_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${CUDA_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_LIBRARIES}")
  endif()
  # MPI
  if (MPI_FOUND AND CHAMELEON_LOOK_FOR_MPI)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    if (MPI_C_LINK_FLAGS)
      if (${MPI_C_LINK_FLAGS} MATCHES "  -")
        string(REGEX REPLACE " -" "-" MPI_C_LINK_FLAGS ${MPI_C_LINK_FLAGS})
      endif()
      list(APPEND REQUIRED_LDFLAGS "${MPI_C_LINK_FLAGS}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_C_LIBRARIES}")
  endif()
  # HWLOC
  if (HWLOC_FOUND)
    if (HWLOC_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${HWLOC_INCLUDE_DIRS}")
    endif()
    if (HWLOC_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${HWLOC_CFLAGS_OTHER}")
    endif()
    if (HWLOC_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${HWLOC_LDFLAGS_OTHER}")
    endif()
    if (HWLOC_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${HWLOC_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${HWLOC_LIBRARIES}")
  endif()
  # TMG
  if (TMG_FOUND)
    if (TMG_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${TMG_INCLUDE_DIRS_DEP}")
    endif()
    if (TMG_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${TMG_CFLAGS_OTHER_DEP}")
    endif()
    if (TMG_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${TMG_LDFLAGS_OTHER_DEP}")
    endif()
    if(TMG_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${TMG_LIBRARY_DIRS_DEP}")
    endif()
    if (TMG_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${TMG_LIBRARIES_DEP}")
    endif()
  endif()
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
  # EXTRA LIBS such that pthread, m, rt
  list(APPEND REQUIRED_LIBS ${CHAMELEON_EXTRA_LIBRARIES})

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
  unset(CHAMELEON_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(MORSE_Init CHAMELEON_WORKS)
  mark_as_advanced(CHAMELEON_WORKS)

  if(CHAMELEON_WORKS)
    # save link with dependencies
    set(CHAMELEON_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(CHAMELEON_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(CHAMELEON_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(CHAMELEON_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(CHAMELEON_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT CHAMELEON_FIND_QUIETLY)
      message(STATUS "Looking for chameleon : test of MORSE_Init fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe CHAMELEON is linked with specific libraries. "
        "Have you tried with COMPONENTS (STARPU/QUARK, CUDA, MPI, FXT)? "
        "See the explanation in FindCHAMELEON.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(CHAMELEON_LIBRARIES)

if (CHAMELEON_LIBRARIES)
  if (CHAMELEON_LIBRARY_DIRS)
    foreach(dir ${CHAMELEON_LIBRARY_DIRS})
      if ("${dir}" MATCHES "chameleon")
        set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET CHAMELEON_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(CHAMELEON_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(CHAMELEON_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of CHAMELEON library" FORCE)
  else()
    set(CHAMELEON_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of CHAMELEON library" FORCE)
  endif()
endif()
mark_as_advanced(CHAMELEON_DIR)
mark_as_advanced(CHAMELEON_DIR_FOUND)

# check that CHAMELEON has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHAMELEON DEFAULT_MSG
  CHAMELEON_LIBRARIES
  CHAMELEON_WORKS)
