###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find MUMPS include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(MUMPS
#               [REQUIRED] # Fail with error if mumps is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  MUMPS depends on the following libraries:
#   - Threads
#   - BLAS
#
#  COMPONENTS are optional libraries MUMPS could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - MPI: to activate detection of the parallel MPI version (default)
#        it looks for Threads, BLAS, MPI and ScaLAPACK libraries
#   - SEQ: to activate detection of sequential version (exclude MPI version)
#        it looks for Threads and BLAS libraries
#   - SCOTCH: to activate detection of MUMPS linked with SCOTCH
#   - PTSCOTCH: to activate detection of MUMPS linked with PTSCOTCH
#   - METIS: to activate detection of MUMPS linked with METIS
#   - PARMETIS: to activate detection of MUMPS linked with PARMETIS
#   - OPENMP: to activate detection of MUMPS linked with OPENMP
#
# This module finds headers and mumps library.
# Results are reported in variables:
#  MUMPS_FOUND             - True if headers and requested libraries were found
#  MUMPS_CFLAGS_OTHER      - mumps compiler flags without headers paths
#  MUMPS_LDFLAGS_OTHER     - mumps linker flags without libraries
#  MUMPS_INCLUDE_DIRS      - mumps include directories
#  MUMPS_LIBRARY_DIRS      - mumps link directories
#  MUMPS_LIBRARIES         - mumps libraries to be linked (absolute path)
#  MUMPS_CFLAGS_OTHER_DEP  - mumps + dependencies compiler flags without headers paths
#  MUMPS_LDFLAGS_OTHER_DEP - mumps + dependencies linker flags without libraries
#  MUMPS_INCLUDE_DIRS_DEP  - mumps + dependencies include directories
#  MUMPS_LIBRARY_DIRS_DEP  - mumps + dependencies link directories
#  MUMPS_LIBRARIES_DEP     - mumps + dependencies libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DMUMPS_DIR=path/to/mumps):
#  MUMPS_DIR              - Where to find the base directory of mumps
# The module can also look for the following environment variables if paths
# are not given as cmake variable: MUMPS_DIR

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
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

if (NOT MUMPS_FOUND)
  set(MUMPS_DIR "" CACHE PATH "Installation directory of MUMPS library")
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "A cache variable, namely MUMPS_DIR, has been set to specify the install directory of MUMPS")
  endif()
endif()

# Set the version to find
set(MUMPS_LOOK_FOR_MPI ON)
set(MUMPS_LOOK_FOR_SEQ OFF)
set(MUMPS_LOOK_FOR_SCOTCH OFF)
set(MUMPS_LOOK_FOR_PTSCOTCH OFF)
set(MUMPS_LOOK_FOR_METIS OFF)
set(MUMPS_LOOK_FOR_PARMETIS OFF)
set(MUMPS_LOOK_FOR_OPENMP OFF)

if( MUMPS_FIND_COMPONENTS )
  foreach( component ${MUMPS_FIND_COMPONENTS} )
    if (${component} STREQUAL "SEQ")
      # means we look for the sequential version of MUMPS (without MPI)
      set(MUMPS_LOOK_FOR_SEQ ON)
      set(MUMPS_LOOK_FOR_MPI OFF)
    endif()
    if (${component} STREQUAL "MPI")
      # means we look for the MPI version of MUMPS (default)
      set(MUMPS_LOOK_FOR_MPI ON)
      set(MUMPS_LOOK_FOR_SEQ OFF)
    endif()
    if (${component} STREQUAL "SCOTCH")
      set(MUMPS_LOOK_FOR_SCOTCH ON)
    endif()
    if (${component} STREQUAL "PTSCOTCH")
      set(MUMPS_LOOK_FOR_PTSCOTCH ON)
    endif()
    if (${component} STREQUAL "METIS")
      set(MUMPS_LOOK_FOR_METIS ON)
    endif()
    if (${component} STREQUAL "PARMETIS")
      set(MUMPS_LOOK_FOR_PARMETIS ON)
    endif()
    if (${component} STREQUAL "OPENMP")
      set(MUMPS_LOOK_FOR_OPENMP ON)
    endif()
  endforeach()
endif()

if (NOT MUMPS_FIND_QUIETLY)
  if (MUMPS_LOOK_FOR_SEQ)
    message(STATUS "Looking for MUMPS - sequential version (without MPI)")
  else()
    message(STATUS "Looking for MUMPS - MPI version -"
      " if you want to force detection of a sequential "
      "version use find_package(MUMPS [REQUIRED] COMPONENTS SEQ [...])")
  endif()
endif()

if (NOT MUMPS_FIND_QUIETLY)
  message(STATUS "Looking for MUMPS - PkgConfig not used")
endif()

# Required dependencies
# ---------------------

if (NOT MUMPS_FIND_QUIETLY)
  message(STATUS "Looking for MUMPS - Try to detect pthread")
endif()
if (MUMPS_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()
set(MUMPS_EXTRA_LIBRARIES "")
if( THREADS_FOUND )
  list(APPEND MUMPS_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# MUMPS depends on BLAS
#----------------------
if (NOT MUMPS_FIND_QUIETLY)
  message(STATUS "Looking for MUMPS - Try to detect BLAS")
endif()
if (MUMPS_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
else()
  find_package(BLAS)
endif()

# Optional dependencies
# ---------------------

# MUMPS may depend on MPI
#------------------------
if (MUMPS_LOOK_FOR_MPI)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect MPI")
  endif()
  # allows to use an external mpi compilation by setting compilers with
  # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
  # at cmake configure
  if(NOT MPI_C_COMPILER)
    set(MPI_C_COMPILER mpicc)
  endif()
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_MPI)
    find_package(MPI REQUIRED)
  else()
    find_package(MPI)
  endif()
  if (MPI_FOUND)
    mark_as_advanced(MPI_LIBRARY)
    mark_as_advanced(MPI_EXTRA_LIBRARY)
  endif()
endif (MUMPS_LOOK_FOR_MPI)

# MUMPS may depend on ScaLAPACK (if MPI version)
#-----------------------------------------------
if (MUMPS_LOOK_FOR_MPI)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect SCALAPACK")
  endif()
  # SCALAPACK is a required dependency if MPI is used
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_MPI)
    find_package(SCALAPACK REQUIRED)
  else()
    find_package(SCALAPACK)
  endif()
endif (MUMPS_LOOK_FOR_MPI)

# MUMPS may depends on SCOTCH
#----------------------------
if (MUMPS_LOOK_FOR_SCOTCH)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect SCOTCH with esmumps")
  endif()
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_SCOTCH)
    find_package(SCOTCH REQUIRED COMPONENTS ESMUMPS)
  else()
    find_package(SCOTCH COMPONENTS ESMUMPS)
  endif()
endif()

# MUMPS may depends on PTSCOTCH
#------------------------------
if (MUMPS_LOOK_FOR_PTSCOTCH)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect PTSCOTCH with esmumps")
  endif()
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_PTSCOTCH)
    find_package(PTSCOTCH REQUIRED COMPONENTS ESMUMPS)
  else()
    find_package(PTSCOTCH COMPONENTS ESMUMPS)
  endif()
endif()

# MUMPS may depends on METIS
#---------------------------
if (MUMPS_LOOK_FOR_METIS)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect METIS")
  endif()
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_METIS)
    find_package(METIS REQUIRED)
  else()
    find_package(METIS)
  endif()
endif()

# MUMPS may depends on PARMETIS
#------------------------------
if (MUMPS_LOOK_FOR_PARMETIS)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect PARMETIS")
  endif()
  if (MUMPS_FIND_REQUIRED AND MUMPS_FIND_REQUIRED_PARMETIS)
    find_package(PARMETIS REQUIRED)
  else()
    find_package(PARMETIS)
  endif()
endif()

# MUMPS may depends on OPENMP
#------------------------------
if (MUMPS_LOOK_FOR_OPENMP)
  if (NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for MUMPS - Try to detect OPENMP")
  endif()
  if (MUMPS_FIND_REQUIRED)
    find_package(OpenMP REQUIRED)
  else()
    find_package(OpenMP)
  endif()
endif()

# Looking for MUMPS
# -----------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_MUMPS_DIR "$ENV{MUMPS_DIR}")
set(ENV_MUMPS_INCDIR "$ENV{MUMPS_INCDIR}")
if(ENV_MUMPS_INCDIR)
  list(APPEND _inc_env "${ENV_MUMPS_INCDIR}")
elseif(ENV_MUMPS_DIR)
  list(APPEND _inc_env "${ENV_MUMPS_DIR}")
  list(APPEND _inc_env "${ENV_MUMPS_DIR}/include")
  list(APPEND _inc_env "${ENV_MUMPS_DIR}/include/mumps")
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

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_MUMPS_LIBDIR "$ENV{MUMPS_LIBDIR}")
if(ENV_MUMPS_LIBDIR)
  list(APPEND _lib_env "${ENV_MUMPS_LIBDIR}")
elseif(ENV_MUMPS_DIR)
  list(APPEND _lib_env "${ENV_MUMPS_DIR}")
  list(APPEND _lib_env "${ENV_MUMPS_DIR}/lib")
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

# Looking for include
# -------------------

# Try to find the mumps header in the given path
# ----------------------------------------------

# create list of headers to find
list(APPEND MUMPS_hdrs_to_find "smumps_c.h;dmumps_c.h;cmumps_c.h;zmumps_c.h")

# call cmake macro to find the header path
if(MUMPS_INCDIR)
  foreach(mumps_hdr ${MUMPS_hdrs_to_find})
    set(MUMPS_${mumps_hdr}_DIRS "MUMPS_${mumps_hdr}_INCLUDE_DIRS-NOTFOUND")
    find_path(MUMPS_${mumps_hdr}_DIRS
      NAMES ${mumps_hdr}
      HINTS ${MUMPS_INCDIR})
  endforeach()
else()
  if(MUMPS_DIR)
    set(MUMPS_${mumps_hdr}_DIRS "MUMPS_${mumps_hdr}_INCLUDE_DIRS-NOTFOUND")
    foreach(mumps_hdr ${MUMPS_hdrs_to_find})
      find_path(MUMPS_${mumps_hdr}_DIRS
        NAMES ${mumps_hdr}
        HINTS ${MUMPS_DIR}
        PATH_SUFFIXES "include")
    endforeach()
  else()
    foreach(mumps_hdr ${MUMPS_hdrs_to_find})
      set(MUMPS_${mumps_hdr}_DIRS "MUMPS_${mumps_hdr}_INCLUDE_DIRS-NOTFOUND")
      find_path(MUMPS_${mumps_hdr}_DIRS
        NAMES ${mumps_hdr}
        HINTS ${_inc_env})
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
# detect which precisions are available
if (MUMPS_smumps_c.h_DIRS)
  mark_as_advanced(MUMPS_smumps_c.h_DIRS)
  set(MUMPS_PREC_S ON)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_smumps_c.h_DIRS}")
else ()
  set(MUMPS_PREC_S OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- smumps_c.h not found")
  endif()
endif()
if (MUMPS_dmumps_c.h_DIRS)
  mark_as_advanced(MUMPS_dmumps_c.h_DIRS)
  set(MUMPS_PREC_D ON)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_dmumps_c.h_DIRS}")
else ()
  set(MUMPS_PREC_D OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- dmumps_c.h not found")
  endif()
endif()
if (MUMPS_cmumps_c.h_DIRS)
  mark_as_advanced(MUMPS_cmumps_c.h_DIRS)
  set(MUMPS_PREC_C ON)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_cmumps_c.h_DIRS}")
else ()
  set(MUMPS_PREC_C OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- cmumps_c.h not found")
  endif()
endif()
if (MUMPS_zmumps_c.h_DIRS)
  mark_as_advanced(MUMPS_zmumps_c.h_DIRS)
  set(MUMPS_PREC_Z ON)
  set(MUMPS_INCLUDE_DIRS "${MUMPS_zmumps_c.h_DIRS}")
else ()
  set(MUMPS_PREC_Z OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- zmumps_c.h not found")
  endif()
endif()


# Looking for lib
# ---------------

# create list of libs to find
set(MUMPS_libs_to_find "mumps_common;pord")
if (MUMPS_LOOK_FOR_SEQ)
  list(APPEND MUMPS_libs_to_find "mpiseq")
endif()
if(MUMPS_PREC_S)
  list(APPEND MUMPS_libs_to_find "smumps")
endif()
if(MUMPS_PREC_D)
  list(APPEND MUMPS_libs_to_find "dmumps")
endif()
if(MUMPS_PREC_C)
  list(APPEND MUMPS_libs_to_find "cmumps")
endif()
if(MUMPS_PREC_Z)
  list(APPEND MUMPS_libs_to_find "zmumps")
endif()

# call cmake macro to find the lib path
if(MUMPS_LIBDIR)
  foreach(mumps_lib ${MUMPS_libs_to_find})
    set(MUMPS_${mumps_lib}_LIBRARY "MUMPS_${mumps_lib}_LIBRARY-NOTFOUND")
    find_library(MUMPS_${mumps_lib}_LIBRARY
      NAMES ${mumps_lib}
      HINTS ${MUMPS_LIBDIR})
  endforeach()
else()
  if(MUMPS_DIR)
    foreach(mumps_lib ${MUMPS_libs_to_find})
      set(MUMPS_${mumps_lib}_LIBRARY "MUMPS_${mumps_lib}_LIBRARY-NOTFOUND")
      find_library(MUMPS_${mumps_lib}_LIBRARY
        NAMES ${mumps_lib}
        HINTS ${MUMPS_DIR}
        PATH_SUFFIXES lib lib32 lib64)
    endforeach()
  else()
    foreach(mumps_lib ${MUMPS_libs_to_find})
      set(MUMPS_${mumps_lib}_LIBRARY "MUMPS_${mumps_lib}_LIBRARY-NOTFOUND")
      find_library(MUMPS_${mumps_lib}_LIBRARY
        NAMES ${mumps_lib}
        HINTS ${_lib_env})
    endforeach()
  endif()
endif()

# If found, add path to cmake variable
# ------------------------------------
set(MUMPS_LIBRARIES "")
set(MUMPS_LIBRARY_DIRS "")
# detect which precisions are available
if (MUMPS_smumps_LIBRARY)
  mark_as_advanced(MUMPS_smumps_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_smumps_LIBRARY}")
  get_filename_component(smumps_lib_path ${MUMPS_smumps_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${smumps_lib_path}")
else ()
  set(MUMPS_PREC_S OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- libsmumps.a not found")
  endif()
endif()
if (MUMPS_dmumps_LIBRARY)
  mark_as_advanced(MUMPS_dmumps_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_dmumps_LIBRARY}")
  get_filename_component(dmumps_lib_path ${MUMPS_dmumps_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${dmumps_lib_path}")
else ()
  set(MUMPS_PREC_D OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- libdmumps.a not found")
  endif()
endif()
if (MUMPS_cmumps_LIBRARY)
  mark_as_advanced(MUMPS_cmumps_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_cmumps_LIBRARY}")
  get_filename_component(cmumps_lib_path ${MUMPS_cmumps_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${cmumps_lib_path}")
else ()
  set(MUMPS_PREC_C OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- libcmumps.a not found")
  endif()
endif()
if (MUMPS_zmumps_LIBRARY)
  mark_as_advanced(MUMPS_zmumps_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_zmumps_LIBRARY}")
  get_filename_component(zmumps_lib_path ${MUMPS_zmumps_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${zmumps_lib_path}")
else ()
  set(MUMPS_PREC_Z OFF)
  if(NOT MUMPS_FIND_QUIETLY)
    message(STATUS "Looking for mumps -- libzmumps.a not found")
  endif()
endif()

# check that one precision arithmetic at least has been discovered
if (NOT MUMPS_PREC_S AND NOT MUMPS_PREC_D AND NOT MUMPS_PREC_C AND NOT MUMPS_PREC_S)
  if (MUMPS_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for mumps -- "
      "no lib[sdcz]mumps.a have been found in ${MUMPS_DIR}/lib when required")
  else()
    if(NOT MUMPS_FIND_QUIETLY)
      message(STATUS "Looking for mumps -- no lib[sdcz]mumps.a have been found")
    endif()
  endif()
endif()
# other MUMPS libraries
if (MUMPS_mumps_common_LIBRARY)
  mark_as_advanced(MUMPS_mumps_common_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_mumps_common_LIBRARY}")
  get_filename_component(mumps_common_lib_path ${MUMPS_mumps_common_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${mumps_common_lib_path}")
else ()
  if (MUMPS_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for mumps -- "
      "libmumps_common.a not found in ${MUMPS_DIR}/lib when required")
  else()
    if(NOT MUMPS_FIND_QUIETLY)
      message(STATUS "Looking for mumps -- libmumps_common.a not found")
    endif()
  endif()
endif()
if (MUMPS_mpiseq_LIBRARY)
  mark_as_advanced(MUMPS_mpiseq_LIBRARY)
  if (MUMPS_LOOK_FOR_SEQ)
    list(APPEND MUMPS_LIBRARIES "${MUMPS_mpiseq_LIBRARY}")
    get_filename_component(mpiseq_lib_path ${MUMPS_mpiseq_LIBRARY} PATH)
    list(APPEND MUMPS_LIBRARY_DIRS "${mpiseq_lib_path}")
    list(APPEND MUMPS_INCLUDE_DIRS "${mpiseq_lib_path}")
  endif()
else ()
  if (MUMPS_FIND_REQUIRED AND MUMPS_LOOK_FOR_SEQ)
    message(FATAL_ERROR "Looking for mumps -- "
      "libmpiseq.a not found in ${MUMPS_DIR}/libseq when required")
  else()
    if(NOT MUMPS_FIND_QUIETLY)
      message(STATUS "Looking for mumps -- libmpiseq.a not found")
    endif()
  endif()
endif()
if (MUMPS_pord_LIBRARY)
  mark_as_advanced(MUMPS_pord_LIBRARY)
  list(APPEND MUMPS_LIBRARIES "${MUMPS_pord_LIBRARY}")
  get_filename_component(pord_lib_path ${MUMPS_pord_LIBRARY} PATH)
  list(APPEND MUMPS_LIBRARY_DIRS "${pord_lib_path}")
else ()
  if (MUMPS_FIND_REQUIRED)
    message(FATAL_ERROR "Looking for mumps -- "
      "libpord.a not found in ${MUMPS_DIR}/lib when required")
  else()
    if(NOT MUMPS_FIND_QUIETLY)
      message(STATUS "Looking for mumps -- libpord.a not found")
    endif()
  endif()
endif()
list(REMOVE_DUPLICATES MUMPS_LIBRARY_DIRS)
list(REMOVE_DUPLICATES MUMPS_INCLUDE_DIRS)

# check a function to validate the find
if(MUMPS_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # MUMPS
  if (MUMPS_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${MUMPS_INCLUDE_DIRS}")
  endif()
  if (MUMPS_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${MUMPS_CFLAGS_OTHER}")
  endif()
  if (MUMPS_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${MUMPS_LDFLAGS_OTHER}")
  endif()
  if (MUMPS_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${MUMPS_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${MUMPS_LIBRARIES}")
  # SCALAPACK
  if (MUMPS_LOOK_FOR_MPI AND SCALAPACK_FOUND)
    if (SCALAPACK_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${SCALAPACK_INCLUDE_DIRS_DEP}")
    endif()
    if (SCALAPACK_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${SCALAPACK_CFLAGS_OTHER_DEP}")
    endif()
    if (SCALAPACK_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${SCALAPACK_LDFLAGS_OTHER_DEP}")
    endif()
    if(SCALAPACK_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${SCALAPACK_LIBRARY_DIRS_DEP}")
    endif()
    if (SCALAPACK_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${SCALAPACK_LIBRARIES_DEP}")
    endif()
  endif()
  # MPI
  if (MUMPS_LOOK_FOR_MPI AND MPI_FOUND)
    if (MPI_C_INCLUDE_PATH)
      list(APPEND REQUIRED_INCDIRS "${MPI_C_INCLUDE_PATH}")
    endif()
    if (MPI_Fortran_LINK_FLAGS)
      if (${MPI_Fortran_LINK_FLAGS} MATCHES "  -")
        string(REGEX REPLACE " -" "-" MPI_Fortran_LINK_FLAGS ${MPI_Fortran_LINK_FLAGS})
      endif()
      list(APPEND REQUIRED_LDFLAGS "${MPI_Fortran_LINK_FLAGS}")
    endif()
    list(APPEND REQUIRED_LIBS "${MPI_Fortran_LIBRARIES}")
  endif()
  # SCOTCH
  if (MUMPS_LOOK_FOR_SCOTCH AND SCOTCH_FOUND)
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
  if (MUMPS_LOOK_FOR_PTSCOTCH AND PTSCOTCH_FOUND)
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
  if (MUMPS_LOOK_FOR_METIS AND METIS_FOUND)
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
  # PARMETIS
  if (MUMPS_LOOK_FOR_PARMETIS AND PARMETIS_FOUND)
    if (PARMETIS_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${PARMETIS_INCLUDE_DIRS_DEP}")
    endif()
    if (PARMETIS_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${PARMETIS_CFLAGS_OTHER_DEP}")
    endif()
    if (PARMETIS_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${PARMETIS_LDFLAGS_OTHER_DEP}")
    endif()
    if(PARMETIS_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${PARMETIS_LIBRARY_DIRS_DEP}")
    endif()
    if (PARMETIS_LIBRARIES_DEP)
      list(APPEND REQUIRED_LIBS "${PARMETIS_LIBRARIES_DEP}")
    endif()
  endif()
  # OpenMP
  if(MUMPS_LOOK_FOR_OPENMP AND OPENMP_FOUND)
    list(APPEND REQUIRED_LDFLAGS "${OpenMP_C_FLAGS}")
  endif()
  # Fortran
  if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    find_library(
      FORTRAN_gfortran_LIBRARY
      NAMES gfortran
      HINTS ${_lib_env}
      )
    mark_as_advanced(FORTRAN_gfortran_LIBRARY)
    if (FORTRAN_gfortran_LIBRARY)
      list(APPEND REQUIRED_LIBS "${FORTRAN_gfortran_LIBRARY}")
    endif()
  elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
    find_library(
      FORTRAN_ifcore_LIBRARY
      NAMES ifcore
      HINTS ${_lib_env}
      )
    mark_as_advanced(FORTRAN_ifcore_LIBRARY)
    if (FORTRAN_ifcore_LIBRARY)
      list(APPEND REQUIRED_LIBS "${FORTRAN_ifcore_LIBRARY}")
    endif()
  endif()
  # EXTRA LIBS such that pthread, m, rt
  list(APPEND REQUIRED_LIBS ${MUMPS_EXTRA_LIBRARIES})

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
  include(CheckFortranFunctionExists)
  unset(MUMPS_PREC_S_WORKS CACHE)
  check_fortran_function_exists(smumps MUMPS_PREC_S_WORKS)
  mark_as_advanced(MUMPS_PREC_S_WORKS)
  unset(MUMPS_PREC_D_WORKS CACHE)
  check_fortran_function_exists(dmumps MUMPS_PREC_D_WORKS)
  mark_as_advanced(MUMPS_PREC_D_WORKS)
  unset(MUMPS_PREC_C_WORKS CACHE)
  check_fortran_function_exists(cmumps MUMPS_PREC_C_WORKS)
  mark_as_advanced(MUMPS_PREC_C_WORKS)
  unset(MUMPS_PREC_Z_WORKS CACHE)
  check_fortran_function_exists(zmumps MUMPS_PREC_Z_WORKS)
  mark_as_advanced(MUMPS_PREC_Z_WORKS)

  set(MUMPS_WORKS FALSE)
  if(MUMPS_PREC_S_WORKS OR MUMPS_PREC_D_WORKS OR MUMPS_PREC_C_WORKS OR MUMPS_PREC_Z_WORKS)
    set(MUMPS_WORKS TRUE)
  endif()

  if(MUMPS_WORKS)
    # save link with dependencies
    set(MUMPS_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(MUMPS_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(MUMPS_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(MUMPS_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(MUMPS_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT MUMPS_FIND_QUIETLY)
      message(STATUS "Looking for MUMPS : test of [sdcz]mumps() fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe MUMPS is linked with specific libraries. "
        "Have you tried with COMPONENTS (MPI/SEQ, SCOTCH, PTSCOTCH, METIS, PARMETIS)? "
        "See the explanation in FindMUMPS.cmake.")
    endif()
  endif()

  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)

endif(MUMPS_LIBRARIES)

if (MUMPS_LIBRARIES)
  list(GET MUMPS_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT MUMPS_LIBRARY_DIRS)
    set(MUMPS_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(MUMPS_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of MUMPS library" FORCE)
  else()
    set(MUMPS_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of MUMPS library" FORCE)
  endif()
endif()
mark_as_advanced(MUMPS_DIR)
mark_as_advanced(MUMPS_DIR_FOUND)

# check that MUMPS has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG
  MUMPS_LIBRARIES
  MUMPS_WORKS)
