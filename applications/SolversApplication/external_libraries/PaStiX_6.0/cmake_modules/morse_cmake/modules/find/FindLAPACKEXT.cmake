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
# - Find LAPACK EXTENDED for MORSE projects: find include dirs and libraries
#
# This module allows to find LAPACK libraries by calling the official FindLAPACK module
# and handles the creation of different library lists whether the user wishes to link
# with a sequential LAPACK or a multihreaded (LAPACK_SEQ_LIBRARIES and LAPACK_PAR_LIBRARIES).
# LAPACK is detected with a FindLAPACK call and if the BLAS vendor is in the following list,
# Intel mkl, ACML then the module tries find the corresponding multithreaded libraries
#
# The following variables have been added to manage links with sequential or multithreaded
# versions:
#  LAPACK_INCLUDE_DIRS    - LAPACK include directories
#  LAPACK_LIBRARY_DIRS    - Link directories for LAPACK libraries
#  LAPACK_SEQ_LIBRARIES   - LAPACK component libraries to be linked (sequential)
#  LAPACK_PAR_LIBRARIES   - LAPACK component libraries to be linked (multithreaded)
#  LAPACKEXT_FOUND        - if a LAPACK has been found
#  LAPACKEXT_LIBRARIES    - Idem LAPACK_LIBRARIES
#  LAPACKEXT_INCLUDE_DIRS - Idem LAPACK_INCLUDE_DIRS
#  LAPACKEXT_LIBRARY_DIRS - Idem LAPACK_LIBRARY_DIRS

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

# Macro to factorize this call. required arguments allows to decide if
# the REQUIRED option must be given to find_package calls
macro(find_package_lapack required)
  if(LAPACKEXT_FIND_REQUIRED AND required)
    if(LAPACKEXT_FIND_QUIETLY)
      find_package(LAPACK REQUIRED QUIET)
    else()
      find_package(LAPACK REQUIRED)
    endif()
  else()
    if(LAPACKEXT_FIND_QUIETLY)
      find_package(LAPACK QUIET)
    else()
      find_package(LAPACK)
    endif()
  endif()
endmacro()

# LAPACKEXT depends on BLAS
if(LAPACKEXT_FIND_QUIETLY)
  find_package(BLASEXT QUIET)
else()
  find_package(BLASEXT)
endif()

if(NOT LAPACKEXT_FIND_QUIETLY)
  message(STATUS "In FindLAPACKEXT")
endif()

if(BLA_VENDOR MATCHES "Intel*")

  ###
  # look for include path if the LAPACK vendor is Intel
  ###

  # gather system include paths
  unset(_inc_env)
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
  list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
  set(ENV_MKLROOT "$ENV{MKLROOT}")
  if (ENV_MKLROOT)
    list(APPEND _inc_env "${ENV_MKLROOT}/include")
  endif()
  list(REMOVE_DUPLICATES _inc_env)

  if (BLAS_DIR)
    set(LAPACK_DIR ${BLAS_DIR})
  endif ()
  if (BLAS_INCDIR)
    set(LAPACK_INCDIR ${BLAS_INCDIR})
  endif ()
  # find mkl.h inside known include paths
  set(LAPACK_mkl_lapack.h_INCLUDE_DIRS "LAPACK_mkl_lapack.h_INCLUDE_DIRS-NOTFOUND")
  if(LAPACK_INCDIR)
    find_path(LAPACK_mkl_lapack.h_INCLUDE_DIRS
      NAMES mkl_lapack.h
      HINTS ${LAPACK_INCDIR})
  else()
    if(LAPACK_DIR)
      find_path(LAPACK_mkl_lapack.h_INCLUDE_DIRS
        NAMES mkl_lapack.h
        HINTS ${LAPACK_DIR}
        PATH_SUFFIXES include)
    else()
      find_path(LAPACK_mkl_lapack.h_INCLUDE_DIRS
        NAMES mkl_lapack.h
        HINTS ${_inc_env})
    endif()
  endif()
  mark_as_advanced(LAPACK_mkl_lapack.h_INCLUDE_DIRS)
  ## Print status if not found
  ## -------------------------
  #if (NOT LAPACK_mkl_lapack.h_INCLUDE_DIRS)
  #    Print_Find_Header_Status(lapack mkl_lapack.h)
  #endif ()
  set(LAPACK_INCLUDE_DIRS "")
  if(LAPACK_mkl_lapack.h_INCLUDE_DIRS)
    list(APPEND LAPACK_INCLUDE_DIRS "${LAPACK_mkl_lapack.h_INCLUDE_DIRS}" )
  endif()

  ###
  # look for libs
  ###

  if (BLA_VENDOR MATCHES "Intel10_64lp*")
    ## look for the sequential version
    set(BLA_VENDOR "Intel10_64lp_seq")
  endif()
  find_package_lapack(0)

  if (LAPACK_FOUND)
    if(BLAS_SEQ_LIBRARIES)
      set(LAPACK_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES}")
    else()
      set(LAPACK_SEQ_LIBRARIES "${LAPACK_SEQ_LIBRARIES-NOTFOUND}")
    endif()
    # if BLAS Intel 10 64 bit -> save sequential and multithreaded versions
    if(BLA_VENDOR MATCHES "Intel10_64lp*")
      if(BLAS_PAR_LIBRARIES)
        set(LAPACK_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES}")
      else()
        set(LAPACK_PAR_LIBRARIES "${LAPACK_PAR_LIBRARIES-NOTFOUND}")
      endif()
    endif()
  endif()

elseif(BLA_VENDOR MATCHES "IBMESSL*")

  ## look for the sequential version
  set(BLA_VENDOR "IBMESSL")
  find_package_lapack(0)

  if (LAPACK_FOUND)
    if(LAPACK_LIBRARIES)
      set(LAPACK_SEQ_LIBRARIES "${LAPACK_LIBRARIES}")
    else()
      set(LAPACK_SEQ_LIBRARIES "${LAPACK_SEQ_LIBRARIES-NOTFOUND}")
    endif()
  endif()

  ## look for the multithreaded version
  set(BLA_VENDOR "IBMESSLMT")
  find_package_lapack(0)

  if (LAPACK_FOUND)
    if(LAPACK_LIBRARIES)
      set(LAPACK_PAR_LIBRARIES "${LAPACK_LIBRARIES}")
    else()
      set(LAPACK_PAR_LIBRARIES "${LAPACK_PAR_LIBRARIES-NOTFOUND}")
    endif()
  endif()

elseif(BLA_VENDOR MATCHES "ACML*")

  ###
  # look for libs
  ###
  find_package_lapack(0)

  if (LAPACK_FOUND)
    if(BLAS_SEQ_LIBRARIES)
      set(LAPACK_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES}")
    else()
      set(LAPACK_SEQ_LIBRARIES "${LAPACK_SEQ_LIBRARIES-NOTFOUND}")
    endif()
    if(BLAS_PAR_LIBRARIES)
      set(LAPACK_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES}")
    else()
      set(LAPACK_PAR_LIBRARIES "${LAPACK_PAR_LIBRARIES-NOTFOUND}")
    endif()
  endif()

else()

  ## look for a sequential version
  # call to the cmake official FindLAPACK module
  # This module sets the following variables:
  #  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
  #    is found
  #  LAPACK_LIBRARIES - uncached list of libraries (using full path name) to
  #    link against to use LAPACK
  #  LAPACK95_LIBRARIES - uncached list of libraries (using full path name)
  #    to link against to use LAPACK95 interface
  #  LAPACK95_FOUND - set to true if a library implementing the LAPACK f95 interface
  #    is found
  #  BLA_STATIC  if set on this determines what kind of linkage we do (static)
  #  BLA_VENDOR  if set checks only the specified vendor, if not set checks
  #     all the possibilities
  #  BLA_F95     if set on tries to find the f95 interfaces for LAPACK/LAPACK
  # Remark: it looks only into paths contained in the system environment variables
  find_package_lapack(0)

  if(LAPACK_FOUND)
    set(LAPACK_SEQ_LIBRARIES "${LAPACK_LIBRARIES}")
  else()
    set(LAPACK_SEQ_LIBRARIES "${LAPACK_SEQ_LIBRARIES-NOTFOUND}")
  endif()
  set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")

endif()

if (LAPACK_SEQ_LIBRARIES)
  set(LAPACK_LIBRARIES "${LAPACK_SEQ_LIBRARIES}")
endif()

# extract libs paths if not given by find_package(LAPACK)
if (NOT LAPACK_LIBRARY_DIRS)
  set(LAPACK_LIBRARY_DIRS "")
  string(REPLACE " " ";" LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
  foreach(lapack_lib ${LAPACK_LIBRARIES})
    if (EXISTS "${lapack_lib}")
      get_filename_component(a_lapack_lib_dir "${lapack_lib}" PATH)
      list(APPEND LAPACK_LIBRARY_DIRS "${a_lapack_lib_dir}" )
    else()
      string(REPLACE "-L" "" lapack_lib "${lapack_lib}")
      if (EXISTS "${lapack_lib}")
        list(APPEND LAPACK_LIBRARY_DIRS "${lapack_lib}" )
      else()
        get_filename_component(a_lapack_lib_dir "${lapack_lib}" PATH)
        if (EXISTS "${a_lapack_lib_dir}")
          list(APPEND LAPACK_LIBRARY_DIRS "${a_lapack_lib_dir}" )
        endif()
      endif()
    endif()
  endforeach()
  if (LAPACK_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES LAPACK_LIBRARY_DIRS)
  endif ()
endif (NOT LAPACK_LIBRARY_DIRS)

# check that LAPACKEXT has been found
# -----------------------------------
include(FindPackageHandleStandardArgs)
if(BLA_VENDOR MATCHES "Intel*")
    if(NOT LAPACKEXT_FIND_QUIETLY)
      message(STATUS "LAPACK found is Intel MKL")
      message(STATUS "LAPACK sequential libraries stored in LAPACK_SEQ_LIBRARIES")
    endif()
    find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
      LAPACK_SEQ_LIBRARIES
      LAPACK_LIBRARY_DIRS
      LAPACK_INCLUDE_DIRS)
    if(BLA_VENDOR MATCHES "Intel10_64lp*" AND LAPACK_PAR_LIBRARIES)
      if(NOT LAPACKEXT_FIND_QUIETLY)
        message(STATUS "LAPACK parallel libraries stored in LAPACK_PAR_LIBRARIES")
      endif()
      find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
        LAPACK_PAR_LIBRARIES)
    endif()
elseif(BLA_VENDOR MATCHES "ACML*")
  if(NOT LAPACKEXT_FIND_QUIETLY)
    message(STATUS "LAPACK found is ACML")
    message(STATUS "LAPACK sequential libraries stored in LAPACK_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
    LAPACK_SEQ_LIBRARIES
    LAPACK_LIBRARY_DIRS)
  if(LAPACK_PAR_LIBRARIES)
    if(NOT LAPACKEXT_FIND_QUIETLY)
      message(STATUS "LAPACK parallel libraries stored in LAPACK_PAR_LIBRARIES")
    endif()
    find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
      LAPACK_PAR_LIBRARIES)
  endif()
elseif(BLA_VENDOR MATCHES "IBMESSL*")
  if(NOT LAPACKEXT_FIND_QUIETLY)
    message(STATUS "LAPACK found is IBMESSL")
    message(STATUS "LAPACK sequential libraries stored in LAPACK_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
    LAPACK_SEQ_LIBRARIES
    LAPACK_LIBRARY_DIRS)
  if(LAPACK_PAR_LIBRARIES)
    if(NOT LAPACKEXT_FIND_QUIETLY)
      message(STATUS "LAPACK parallel libraries stored in LAPACK_PAR_LIBRARIES")
    endif()
    find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
      LAPACK_PAR_LIBRARIES)
  endif()
else()
  if(NOT LAPACKEXT_FIND_QUIETLY)
    message(STATUS "LAPACK sequential libraries stored in LAPACK_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(LAPACKEXT DEFAULT_MSG
    LAPACK_SEQ_LIBRARIES
    LAPACK_LIBRARY_DIRS)
endif()

if (LAPACK_LIBRARIES)
  set(LAPACKEXT_LIBRARIES ${LAPACK_LIBRARIES})
endif()
if (LAPACK_INCLUDE_DIRS)
  set(LAPACKEXT_INCLUDE_DIRS ${LAPACK_INCLUDE_DIRS})
endif()
if (LAPACK_LIBRARY_DIRS)
  set(LAPACKEXT_LIBRARY_DIRS ${LAPACK_LIBRARY_DIRS})
endif()
