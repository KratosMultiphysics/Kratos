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
# - Find PASTIX include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PASTIX
#               [REQUIRED] # Fail with error if pastix is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  PASTIX depends on the following libraries:
#   - Threads, m, rt
#   - MPI
#   - HWLOC
#   - BLAS
#   - SPM
#
#  COMPONENTS are optional libraries PASTIX could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - MPI: to activate detection of the parallel MPI version (default)
#        it looks for Threads, HWLOC, BLAS, MPI and ScaLAPACK libraries
#   - SEQ: to activate detection of the sequential version (exclude MPI version)
#   - PARSEC: to activate detection of PASTIX with PARSEC
#   - STARPU: to activate detection of StarPU version
#   it looks for MPI version of StarPU (default behaviour)
#   if SEQ and STARPU are given, it looks for a StarPU without MPI
#   - STARPU_CUDA: to activate detection of StarPU with CUDA
#   - STARPU_FXT: to activate detection of StarPU with FxT
#   - SCOTCH: to activate detection of PASTIX linked with SCOTCH
#   - PTSCOTCH: to activate detection of PASTIX linked with SCOTCH
#   - METIS: to activate detection of PASTIX linked with SCOTCH
#   - FORTRAN: to provide Fortran libraries of PASTIX in LIBRARIES
#
# This module finds headers and pastix library.
# Results are reported in variables:
#  PASTIX_FOUND             - True if headers and requested libraries were found
#  PASTIX_CFLAGS_OTHER      - pastix compiler flags without headers paths
#  PASTIX_LDFLAGS_OTHER     - pastix linker flags without libraries
#  PASTIX_INCLUDE_DIRS      - pastix include directories
#  PASTIX_LIBRARY_DIRS      - pastix link directories
#  PASTIX_LIBRARIES         - pastix libraries to be linked (absolute path)
#  PASTIX_LIBRARIES_FORTRAN - pastix Fortran libraries to be linked (absolute path)
#  PASTIX_CFLAGS_OTHER_DEP  - pastix + dependencies compiler flags without headers paths
#  PASTIX_LDFLAGS_OTHER_DEP - pastix + dependencies linker flags without libraries
#  PASTIX_INCLUDE_DIRS_DEP  - pastix + dependencies include directories
#  PASTIX_LIBRARY_DIRS_DEP  - pastix + dependencies link directories
#  PASTIX_LIBRARIES_DEP     - pastix + dependencies libraries
#  PASTIX_INTSIZE           - Number of octets occupied by a pastix_int_t
#  PASTIX_VERSION_MAJOR     - pastix Version, first (major) number
#
#  PASTIX_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = PASTIX
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
# options at configure (ex: cmake path/to/project -DPASTIX_DIR=path/to/pastix):
#  PASTIX_DIR              - Where to find the base directory of pastix
#  PASTIX_INCDIR           - Where to find the header files
#  PASTIX_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PASTIX_DIR, PASTIX_INCDIR, PASTIX_LIBDIR

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

if (NOT PASTIX_FOUND)
  set(PASTIX_DIR "" CACHE PATH "Installation directory of PASTIX library")
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PASTIX_DIR, has been set to specify the install directory of PASTIX")
  endif()
endif()

# Set the version to find
set(PASTIX_LOOK_FOR_MPI ON)
set(PASTIX_LOOK_FOR_SEQ OFF)
set(PASTIX_LOOK_FOR_STARPU OFF)
set(PASTIX_LOOK_FOR_PARSEC OFF)
set(PASTIX_LOOK_FOR_STARPU_CUDA OFF)
set(PASTIX_LOOK_FOR_STARPU_FXT OFF)
set(PASTIX_LOOK_FOR_SCOTCH ON)
set(PASTIX_LOOK_FOR_PTSCOTCH OFF)
set(PASTIX_LOOK_FOR_METIS OFF)
set(PASTIX_LOOK_FOR_FORTRAN OFF)

if( PASTIX_FIND_COMPONENTS )
  foreach( component ${PASTIX_FIND_COMPONENTS} )
    if (${component} STREQUAL "SEQ")
      # means we look for the sequential version of PaStiX (without MPI)
      set(PASTIX_LOOK_FOR_SEQ ON)
      set(PASTIX_LOOK_FOR_MPI OFF)
    endif()
    if (${component} STREQUAL "MPI")
      # means we look for the MPI version of PaStiX (default)
      set(PASTIX_LOOK_FOR_SEQ OFF)
      set(PASTIX_LOOK_FOR_MPI ON)
    endif()
    if (${component} STREQUAL "PARSEC")
      set(PASTIX_LOOK_FOR_PARSEC ON)
    endif()
    if (${component} STREQUAL "STARPU")
      # means we look for PaStiX with StarPU
      set(PASTIX_LOOK_FOR_STARPU ON)
    endif()
    if (${component} STREQUAL "STARPU_CUDA")
      # means we look for PaStiX with StarPU + CUDA
      set(PASTIX_LOOK_FOR_STARPU ON)
      set(PASTIX_LOOK_FOR_STARPU_CUDA ON)
    endif()
    if (${component} STREQUAL "STARPU_FXT")
      # means we look for PaStiX with StarPU + FxT
      set(PASTIX_LOOK_FOR_STARPU_FXT ON)
    endif()
    if (${component} STREQUAL "SCOTCH")
      set(PASTIX_LOOK_FOR_SCOTCH ON)
    endif()
    if (${component} STREQUAL "SCOTCH")
      set(PASTIX_LOOK_FOR_PTSCOTCH ON)
    endif()
    if (${component} STREQUAL "METIS")
      set(PASTIX_LOOK_FOR_METIS ON)
    endif()
    if (${component} STREQUAL "FORTRAN")
      set(PASTIX_LOOK_FOR_FORTRAN ON)
    endif()
  endforeach()
endif()

# Dependencies detection
# ----------------------

# Required dependencies
# ---------------------

if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect pthread")
endif()
if (PASTIX_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()
set(PASTIX_EXTRA_LIBRARIES "")
if( THREADS_FOUND )
  list(APPEND PASTIX_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect libm")
endif()
set(PASTIX_M_LIBRARIES "")
if(UNIX OR WIN32)
  find_library(
    PASTIX_M_m_LIBRARY
    NAMES m
    )
  mark_as_advanced(PASTIX_M_m_LIBRARY)
  if (PASTIX_M_m_LIBRARY)
    list(APPEND PASTIX_M_LIBRARIES "${PASTIX_M_m_LIBRARY}")
    list(APPEND PASTIX_EXTRA_LIBRARIES "${PASTIX_M_m_LIBRARY}")
  else()
    if (PASTIX_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find libm on your system."
        "Are you sure to a have a C compiler installed?")
    endif()
  endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect librt")
endif()
set(PASTIX_RT_LIBRARIES "")
if(UNIX AND NOT APPLE)
  find_library(
    PASTIX_RT_rt_LIBRARY
    NAMES rt
    )
  mark_as_advanced(PASTIX_RT_rt_LIBRARY)
  if (PASTIX_RT_rt_LIBRARY)
    list(APPEND PASTIX_RT_LIBRARIES "${PASTIX_RT_rt_LIBRARY}")
    list(APPEND PASTIX_EXTRA_LIBRARIES "${PASTIX_RT_rt_LIBRARY}")
  else()
    if (PASTIX_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
  endif()
endif()

# PASTIX depends on HWLOC
#------------------------
if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect HWLOC")
endif()
if (PASTIX_FIND_REQUIRED)
  find_package(HWLOC REQUIRED)
else()
  find_package(HWLOC)
endif()

# PASTIX depends on BLAS
#-----------------------
if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect BLAS")
endif()
if (PASTIX_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
else()
  find_package(BLAS)
endif()

# PASTIX depends on SPM
#------------------------
if (NOT PASTIX_FIND_QUIETLY)
  message(STATUS "Looking for PASTIX - Try to detect SPM")
endif()
set(SPM_COMPONENTS "")
if(PASTIX_LOOK_FOR_FORTRAN)
  list(APPEND SPM_COMPONENTS "FORTRAN" )
endif()
if (PASTIX_FIND_REQUIRED)
  find_package(SPM REQUIRED COMPONENTS ${SPM_COMPONENTS})
else()
  find_package(SPM COMPONENTS ${SPM_COMPONENTS})
endif()

# Optional dependencies
# ---------------------

# PASTIX may depend on MPI
#-------------------------
if (PASTIX_LOOK_FOR_MPI)
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect MPI")
  endif()
  # allows to use an external mpi compilation by setting compilers with
  # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
  # at cmake configure
  if(NOT MPI_C_COMPILER)
    set(MPI_C_COMPILER mpicc)
  endif()
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_MPI)
    find_package(MPI REQUIRED)
  else()
    find_package(MPI)
  endif()
  if (MPI_FOUND)
    mark_as_advanced(MPI_LIBRARY)
    mark_as_advanced(MPI_EXTRA_LIBRARY)
  endif()
endif (PASTIX_LOOK_FOR_MPI)

# PASTIX may depends on SCOTCH
#-----------------------------
if (PASTIX_LOOK_FOR_PARSEC)
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect PARSEC")
  endif()
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_PARSEC)
    find_package(PARSEC REQUIRED)
  else()
    find_package(PARSEC)
  endif()
endif()

# PASTIX may depend on STARPU
#----------------------------
if(PASTIX_LOOK_FOR_STARPU)

  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect StarPU")
  endif()

  set(PASTIX_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")

  # create list of components in order to make a single call to find_package(starpu...)
  # we explicitly need a StarPU version built with hwloc
  set(STARPU_COMPONENT_LIST "HWLOC")

  # StarPU may depend on MPI
  # allows to use an external mpi compilation by setting compilers with
  # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
  # at cmake configure
  if (PASTIX_LOOK_FOR_MPI)
    if(NOT MPI_C_COMPILER)
      set(MPI_C_COMPILER mpicc)
    endif()
    list(APPEND STARPU_COMPONENT_LIST "MPI")
  endif()
  if (PASTIX_LOOK_FOR_STARPU_CUDA)
    list(APPEND STARPU_COMPONENT_LIST "CUDA")
  endif()
  if (PASTIX_LOOK_FOR_STARPU_FXT)
    list(APPEND STARPU_COMPONENT_LIST "FXT")
  endif()
  # set the list of optional dependencies we may discover
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_STARPU)
    find_package(STARPU ${PASTIX_STARPU_VERSION} REQUIRED
      COMPONENTS ${STARPU_COMPONENT_LIST})
  else()
    find_package(STARPU ${PASTIX_STARPU_VERSION}
      COMPONENTS ${STARPU_COMPONENT_LIST})
  endif()

endif(PASTIX_LOOK_FOR_STARPU)

# PASTIX may depends on SCOTCH
#-----------------------------
if (PASTIX_LOOK_FOR_SCOTCH)
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect SCOTCH")
  endif()
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_SCOTCH)
    find_package(SCOTCH REQUIRED)
  else()
    find_package(SCOTCH)
  endif()
endif()

# PASTIX may depends on PTSCOTCH
#-------------------------------
if (PASTIX_LOOK_FOR_PTSCOTCH)
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect PTSCOTCH")
  endif()
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_PTSCOTCH)
    find_package(PTSCOTCH REQUIRED)
  else()
    find_package(PTSCOTCH)
  endif()
endif()

# PASTIX may depends on METIS
#----------------------------
if (PASTIX_LOOK_FOR_METIS)
  if (NOT PASTIX_FIND_QUIETLY)
    message(STATUS "Looking for PASTIX - Try to detect METIS")
  endif()
  if (PASTIX_FIND_REQUIRED AND PASTIX_FIND_REQUIRED_METIS)
    find_package(METIS REQUIRED)
  else()
    find_package(METIS)
  endif()
endif()

set(ENV_PASTIX_DIR "$ENV{PASTIX_DIR}")
set(ENV_PASTIX_INCDIR "$ENV{PASTIX_INCDIR}")
set(ENV_PASTIX_LIBDIR "$ENV{PASTIX_LIBDIR}")
set(PASTIX_GIVEN_BY_USER "FALSE")
if ( PASTIX_DIR OR ( PASTIX_INCDIR AND PASTIX_LIBDIR) OR ENV_PASTIX_DIR OR (ENV_PASTIX_INCDIR AND ENV_PASTIX_LIBDIR) )
  set(PASTIX_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if(PKG_CONFIG_EXECUTABLE AND NOT PASTIX_GIVEN_BY_USER)

  pkg_search_module(PASTIX pastix)

  if (NOT PASTIX_FIND_QUIETLY)
    if (PASTIX_FOUND AND PASTIX_LIBRARIES)
      message(STATUS "Looking for PASTIX - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for PASTIX - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing pastix.pc"
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  if (PASTIX_FOUND AND PASTIX_LIBRARIES)
    set(PASTIX_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(PASTIX)
  else()
    set(PASTIX_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

  set(PASTIX_VERSION_MAJOR)
  foreach(_lib ${PASTIX_LIBRARIES})
    if (_lib MATCHES "pastix_murge")
      set(PASTIX_VERSION_MAJOR "5")
    endif()
    if (_lib MATCHES "pastix_kernels")
      set(PASTIX_VERSION_MAJOR "6")
    endif()
  endforeach()

endif(PKG_CONFIG_EXECUTABLE AND NOT PASTIX_GIVEN_BY_USER)

# Error if pastix required and no partitioning lib found
if (PASTIX_FIND_REQUIRED AND NOT SCOTCH_FOUND AND NOT PTSCOTCH_FOUND AND NOT METIS_FOUND)
  message(FATAL_ERROR "Could NOT find any partitioning library on your system (install scotch, ptscotch or metis)")
endif()

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PASTIX_FOUND) OR (PASTIX_GIVEN_BY_USER) )

  # Looking for PaStiX
  # ------------------

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_PASTIX_DIR "$ENV{PASTIX_DIR}")
  set(ENV_PASTIX_INCDIR "$ENV{PASTIX_INCDIR}")
  if(ENV_PASTIX_INCDIR)
    list(APPEND _inc_env "${ENV_PASTIX_INCDIR}")
  elseif(ENV_PASTIX_DIR)
    list(APPEND _inc_env "${ENV_PASTIX_DIR}")
    list(APPEND _inc_env "${ENV_PASTIX_DIR}/include")
    list(APPEND _inc_env "${ENV_PASTIX_DIR}/include/pastix")
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


  # Try to find the pastix header in the given paths
  # ---------------------------------------------------
  # call cmake macro to find the header path
  if(PASTIX_INCDIR)
    set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
    find_path(PASTIX_pastix.h_DIRS
      NAMES pastix.h
      HINTS ${PASTIX_INCDIR})
  else()
    if(PASTIX_DIR)
      set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
      find_path(PASTIX_pastix.h_DIRS
        NAMES pastix.h
        HINTS ${PASTIX_DIR}
        PATH_SUFFIXES "include" "include/pastix")
    else()
      set(PASTIX_pastix.h_DIRS "PASTIX_pastix.h_DIRS-NOTFOUND")
      find_path(PASTIX_pastix.h_DIRS
        NAMES pastix.h
        HINTS ${_inc_env}
        PATH_SUFFIXES "pastix")
    endif()
  endif()
  mark_as_advanced(PASTIX_pastix.h_DIRS)

  # If found, add path to cmake variable
  # ------------------------------------
  if (PASTIX_pastix.h_DIRS)
    set(PASTIX_INCLUDE_DIRS "${PASTIX_pastix.h_DIRS}")
  else ()
    set(PASTIX_INCLUDE_DIRS "PASTIX_INCLUDE_DIRS-NOTFOUND")
    if(NOT PASTIX_FIND_QUIETLY)
      message(STATUS "Looking for pastix -- pastix.h not found")
    endif()
  endif()


  # Looking for lib
  # ---------------

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_PASTIX_LIBDIR "$ENV{PASTIX_LIBDIR}")
  if(ENV_PASTIX_LIBDIR)
    list(APPEND _lib_env "${ENV_PASTIX_LIBDIR}")
  elseif(ENV_PASTIX_DIR)
    list(APPEND _lib_env "${ENV_PASTIX_DIR}")
    list(APPEND _lib_env "${ENV_PASTIX_DIR}/lib")
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

  # Try to find the pastix lib in the given paths
  # ------------------------------------------------

  # create list of libs to find
  set(PASTIX_libs_to_find_v5 "pastix_murge;pastix")
  set(PASTIX_libs_to_find_v6 "pastix;pastix_kernels")
  if(PASTIX_LOOK_FOR_PARSEC)
    list(INSERT 1 PASTIX_libs_to_find_v6 "pastix_parsec")
  endif()
  if(PASTIX_LOOK_FOR_STARPU)
    list(INSERT 1 PASTIX_libs_to_find_v6 "pastix_starpu")
  endif()
  foreach(v 5 6)
    set(V_FOUND TRUE)
    set(PASTIX_libs_to_find ${PASTIX_libs_to_find_v${v}})
    message(STATUS "Looking for pastix -- Testing is version ${v}")
    # call cmake macro to find the lib path
    foreach(pastix_lib ${PASTIX_libs_to_find})
      set(PASTIX_${pastix_lib}_LIBRARY "PASTIX_${pastix_lib}_LIBRARY-NOTFOUND")
      if(PASTIX_LIBDIR)
        find_library(PASTIX_${pastix_lib}_LIBRARY
          NAMES ${pastix_lib}
          HINTS ${PASTIX_LIBDIR})
      else()
        if(PASTIX_DIR)
          find_library(PASTIX_${pastix_lib}_LIBRARY
            NAMES ${pastix_lib}
            HINTS ${PASTIX_DIR}
            PATH_SUFFIXES lib lib32 lib64)
        else()
          find_library(PASTIX_${pastix_lib}_LIBRARY
            NAMES ${pastix_lib}
            HINTS ${_lib_env})
        endif()
      endif()
      if (NOT PASTIX_${pastix_lib}_LIBRARY)
        set (V_FOUND FALSE)
        break()
      endif()
    endforeach()
    if (V_FOUND)
      set(PASTIX_VERSION_MAJOR ${v})
      message(STATUS "Looking for pastix -- major version detected: ${PASTIX_VERSION_MAJOR}")
      if (PASTIX_FIND_VERSION)
        if (NOT ${PASTIX_VERSION_MAJOR} EQUAL ${PASTIX_FIND_VERSION_MAJOR})
          message(FATAL_ERROR "Looking for pastix -- PaStiX version found differ from required")
        endif()
      endif()
      break()
    endif()
  endforeach()

  # If found, add path to cmake variable
  # ------------------------------------
  foreach(pastix_lib ${PASTIX_libs_to_find})

    get_filename_component(${pastix_lib}_lib_path ${PASTIX_${pastix_lib}_LIBRARY} PATH)
    # set cmake variables (respects naming convention)
    if (PASTIX_LIBRARIES)
      list(APPEND PASTIX_LIBRARIES "${PASTIX_${pastix_lib}_LIBRARY}")
    else()
      set(PASTIX_LIBRARIES "${PASTIX_${pastix_lib}_LIBRARY}")
    endif()
    if (PASTIX_LIBRARY_DIRS)
      list(APPEND PASTIX_LIBRARY_DIRS "${${pastix_lib}_lib_path}")
    else()
      set(PASTIX_LIBRARY_DIRS "${${pastix_lib}_lib_path}")
    endif()
    mark_as_advanced(PASTIX_${pastix_lib}_LIBRARY)

  endforeach(pastix_lib ${PASTIX_libs_to_find})

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PASTIX_FOUND) OR (PASTIX_GIVEN_BY_USER) )


# check a function to validate the find
if(PASTIX_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PASTIX
  if (PASTIX_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${PASTIX_INCLUDE_DIRS}")
  endif()
  if (PASTIX_CFLAGS_OTHER)
    list(APPEND REQUIRED_FLAGS "${PASTIX_CFLAGS_OTHER}")
  endif()
  if (PASTIX_LDFLAGS_OTHER)
    list(APPEND REQUIRED_LDFLAGS "${PASTIX_LDFLAGS_OTHER}")
  endif()
  foreach(libdir ${PASTIX_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${PASTIX_LIBRARIES}")
  # PARSEC
  if (PASTIX_LOOK_FOR_PARSEC AND PARSEC_FOUND)
    if (PARSEC_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${PARSEC_INCLUDE_DIRS_DEP}")
    endif()
    if (PARSEC_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${PARSEC_CFLAGS_OTHER_DEP}")
    endif()
    if (PARSEC_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${PARSEC_LDFLAGS_OTHER_DEP}")
    endif()
    if(PARSEC_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${PARSEC_LIBRARY_DIRS_DEP}")
    endif()
    if (PARSEC_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${PARSEC_LIBRARIES_DEP}")
    endif()
  endif()
  # STARPU
  if (PASTIX_LOOK_FOR_STARPU AND STARPU_FOUND)
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
  # CUDA
  if (PASTIX_LOOK_FOR_STARPU_CUDA AND CUDA_FOUND)
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
  if (PASTIX_LOOK_FOR_MPI AND MPI_FOUND)
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
  # BLAS
  if (BLAS_FOUND)
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
  endif()
  # SPM
  if (SPM_FOUND)
    if (SPM_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${SPM_INCLUDE_DIRS}")
    endif()
    if (SPM_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${SPM_CFLAGS_OTHER}")
    endif()
    if (SPM_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${SPM_LDFLAGS_OTHER}")
    endif()
    if (SPM_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${SPM_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${SPM_LIBRARIES}")
  endif()
  # SCOTCH
  if (PASTIX_LOOK_FOR_SCOTCH AND SCOTCH_FOUND)
    if (SCOTCH_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${SCOTCH_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${SCOTCH_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${SCOTCH_LIBRARIES}")
  endif()
  # PTSCOTCH
  if (PASTIX_LOOK_FOR_PTSCOTCH AND PTSCOTCH_FOUND)
    if (PTSCOTCH_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${PTSCOTCH_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${PTSCOTCH_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${PTSCOTCH_LIBRARIES}")
  endif()
  # METIS
  if (PASTIX_LOOK_FOR_METIS AND METIS_FOUND)
    if (METIS_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${METIS_INCLUDE_DIRS}")
    endif()
    foreach(libdir ${METIS_LIBRARY_DIRS})
      if (libdir)
        list(APPEND REQUIRED_LIBDIRS "${libdir}")
      endif()
    endforeach()
    list(APPEND REQUIRED_LIBS "${METIS_LIBRARIES}")
  endif()
  # EXTRA LIBS such that pthread, m, rt
  list(APPEND REQUIRED_LIBS ${PASTIX_EXTRA_LIBRARIES})

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
  unset(PASTIX_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(pastix PASTIX_WORKS)
  mark_as_advanced(PASTIX_WORKS)

  if(PASTIX_WORKS)
    # save link with dependencies
    set(PASTIX_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(PASTIX_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(PASTIX_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(PASTIX_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(PASTIX_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT PASTIX_FIND_QUIETLY)
      message(STATUS "Looking for PASTIX : test of pastix() fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe PASTIX is linked with specific libraries. "
        "Have you tried with COMPONENTS (MPI/SEQ, PARSEC, STARPU, STARPU_CUDA, SCOTCH, PTSCOTCH, METIS)? "
        "See the explanation in FindPASTIX.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(PASTIX_LIBRARIES)


# Check the size of pastix_int_t
# ---------------------------------
set(CMAKE_REQUIRED_INCLUDES ${PASTIX_INCLUDE_DIRS_DEP})

include(CheckCSourceRuns)
#mpi.h should be included by pastix.h directly
if ((PASTIX_LOOK_FOR_MPI AND MPI_FOUND) AND (${PASTIX_VERSION_MAJOR} EQUAL 5))
  set(PASTIX_C_TEST_PASTIX_INT_4 "
#include <mpi.h>
#include <pastix.h>
int main(int argc, char **argv) {
  if (sizeof(pastix_int_t) == 4)
    return 0;
  else
    return 1;
}
")

  set(PASTIX_C_TEST_PASTIX_INT_8 "
#include <mpi.h>
#include <pastix.h>
int main(int argc, char **argv) {
  if (sizeof(pastix_int_t) == 8)
    return 0;
  else
    return 1;
}
")
else()
  set(PASTIX_C_TEST_PASTIX_INT_4 "
#include <pastix.h>
int main(int argc, char **argv) {
  if (sizeof(pastix_int_t) == 4)
    return 0;
  else
    return 1;
}
")

  set(PASTIX_C_TEST_PASTIX_INT_8 "
#include <pastix.h>
int main(int argc, char **argv) {
  if (sizeof(pastix_int_t) == 8)
    return 0;
  else
    return 1;
}
")
endif()

unset(PASTIX_INT_4 CACHE)
unset(PASTIX_INT_8 CACHE)
check_c_source_runs("${PASTIX_C_TEST_PASTIX_INT_4}" PASTIX_INT_4)
check_c_source_runs("${PASTIX_C_TEST_PASTIX_INT_8}" PASTIX_INT_8)
if(NOT PASTIX_INT_4)
  if(NOT PASTIX_INT_8)
    set(PASTIX_INTSIZE -1)
  else()
    set(PASTIX_INTSIZE 8)
  endif()
else()
  set(PASTIX_INTSIZE 4)
endif()
set(CMAKE_REQUIRED_INCLUDES "")


if (PASTIX_LIBRARIES)
  list(GET PASTIX_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT PASTIX_LIBRARY_DIRS)
    set(PASTIX_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PASTIX_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PASTIX library" FORCE)
  else()
    set(PASTIX_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PASTIX library" FORCE)
  endif()
endif()
mark_as_advanced(PASTIX_DIR)
mark_as_advanced(PASTIX_DIR_FOUND)

if (PASTIX_LOOK_FOR_FORTRAN)
  find_library(PASTIX_pastixf_LIBRARY
               NAMES pastixf
               HINTS ${PASTIX_LIBRARY_DIRS})
  mark_as_advanced(PASTIX_pastixf_LIBRARY)
  if (PASTIX_pastixf_LIBRARY)
    set(PASTIX_LIBRARIES_FORTRAN ${PASTIX_pastixf_LIBRARY})
    list(INSERT PASTIX_LIBRARIES 0 "${PASTIX_LIBRARIES_FORTRAN}")
    list(INSERT PASTIX_LIBRARIES_DEP 0 "${PASTIX_LIBRARIES_FORTRAN}")
  else()
    set(PASTIX_LIBRARIES_FORTRAN "PASTIX_LIBRARIES_FORTRAN-NOTFOUND")
  endif()
endif()

# check that PASTIX has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if (PASTIX_LOOK_FOR_FORTRAN)
  find_package_handle_standard_args(PASTIX DEFAULT_MSG
    PASTIX_LIBRARIES_FORTRAN
    PASTIX_LIBRARIES
    PASTIX_WORKS)
else()
  find_package_handle_standard_args(PASTIX DEFAULT_MSG
    PASTIX_LIBRARIES
    PASTIX_WORKS)
endif()
