###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2017 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find BLAS EXTENDED for MORSE projects: find include dirs and libraries
#
# This module allows to find BLAS libraries by calling the official FindBLAS module
# and handles the creation of different library lists whether the user wishes to link
# with a sequential BLAS or a multihreaded (BLAS_SEQ_LIBRARIES and BLAS_PAR_LIBRARIES).
# BLAS is detected with a FindBLAS call then if the BLAS vendor is Intel10_64lp, ACML
# or IBMESSLMT then the module attempts to find the corresponding multithreaded libraries.
#
# The following variables have been added to manage links with sequential or multithreaded
# versions:
#  BLAS_INCLUDE_DIRS      - BLAS include directories
#  BLAS_LIBRARY_DIRS      - Link directories for BLAS libraries
#  BLAS_SEQ_LIBRARIES     - BLAS component libraries to be linked (sequential)
#  BLAS_SEQ_CFLAGS_OTHER  - compiler flags without headers paths
#  BLAS_SEQ_LDFLAGS_OTHER - linker flags without libraries
#  BLAS_PAR_LIBRARIES     - BLAS component libraries to be linked (multithreaded)
#  BLAS_PAR_CFLAGS_OTHER  - compiler flags without headers paths
#  BLAS_PAR_LDFLAGS_OTHER - linker flags without libraries
#  BLASEXT_FOUND          - if a BLAS has been found
#  BLASEXT_LIBRARIES      - Idem BLAS_LIBRARIES
#  BLASEXT_INCLUDE_DIRS   - Idem BLAS_INCLUDE_DIRS
#  BLASEXT_LIBRARY_DIRS   - Idem BLAS_LIBRARY_DIRS

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013-2017 Florent Pruvost
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
macro(find_package_blas required)
  if(BLASEXT_FIND_REQUIRED AND required)
    if(BLASEXT_FIND_QUIETLY)
      find_package(BLAS REQUIRED QUIET)
    else()
      find_package(BLAS REQUIRED)
    endif()
  else()
    if(BLASEXT_FIND_QUIETLY)
      find_package(BLAS QUIET)
    else()
      find_package(BLAS)
    endif()
  endif()
endmacro()

# add a cache variable to let the user specify the BLAS vendor
set(BLA_VENDOR "" CACHE STRING "list of possible BLAS vendor:
    Open, Eigen, Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, IBMESSLMT,
    Intel10_32 (intel mkl v10 32 bit),
    Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model),
    Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),
    Intel( older versions of mkl 32 and 64 bit),
    ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")

if(NOT BLASEXT_FIND_QUIETLY)
  message(STATUS "In FindBLASEXT")
  message(STATUS "If you want to force the use of one specific library, "
    "\n   please specify the BLAS vendor by setting -DBLA_VENDOR=blas_vendor_name"
    "\n   at cmake configure.")
  message(STATUS "List of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, "
    "\n   DXML, SunPerf, SCSL, SGIMATH, IBMESSL, IBMESSLMT, Intel10_32 (intel mkl v10 32 bit),"
    "\n   Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model),"
    "\n   Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
    "\n   Intel( older versions of mkl 32 and 64 bit),"
    "\n   ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
endif()

# First blas detection in order to decide if we should look for a
# multitheaded version
find_package_blas(0)

# detect the cases where SEQ and PAR libs are handled
if(BLA_VENDOR STREQUAL "All" AND
    (BLAS_mkl_core_LIBRARY OR BLAS_mkl_core_dll_LIBRARY)
    )
  set(BLA_VENDOR "Intel")
  if(BLAS_mkl_intel_LIBRARY)
    set(BLA_VENDOR "Intel10_32")
  endif()
  if(BLAS_mkl_intel_lp64_LIBRARY)
    set(BLA_VENDOR "Intel10_64lp")
  endif()
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "A BLAS library has been found (${BLAS_LIBRARIES}) but we"
      "\n   have also potentially detected some multithreaded BLAS libraries from the MKL."
      "\n   We try to find both libraries lists (Sequential/Multithreaded).")
  endif()
  set(BLAS_FOUND "")
elseif(BLA_VENDOR STREQUAL "All" AND BLAS_acml_LIBRARY)
  set(BLA_VENDOR "ACML")
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "A BLAS library has been found (${BLAS_LIBRARIES}) but we"
      "\n   have also potentially detected some multithreaded BLAS libraries from the ACML."
      "\n   We try to find both libraries lists (Sequential/Multithreaded).")
  endif()
  set(BLAS_FOUND "")
elseif(BLA_VENDOR STREQUAL "All" AND BLAS_essl_LIBRARY)
  set(BLA_VENDOR "IBMESSL")
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "A BLAS library has been found (${BLAS_LIBRARIES}) but we"
      "\n   have also potentially detected some multithreaded BLAS libraries from the ESSL."
      "\n   We try to find both libraries lists (Sequential/Multithreaded).")
  endif()
  set(BLAS_FOUND "")
endif()

# Intel case
if(BLA_VENDOR MATCHES "Intel*")

  ###
  # look for include path if the BLAS vendor is Intel
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

  # find mkl.h inside known include paths
  set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
  if(BLAS_INCDIR)
    set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
    find_path(BLAS_mkl.h_INCLUDE_DIRS
      NAMES mkl.h
      HINTS ${BLAS_INCDIR})
  else()
    if(BLAS_DIR)
      set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
      find_path(BLAS_mkl.h_INCLUDE_DIRS
        NAMES mkl.h
        HINTS ${BLAS_DIR}
        PATH_SUFFIXES include)
    else()
      set(BLAS_mkl.h_INCLUDE_DIRS "BLAS_mkl.h_INCLUDE_DIRS-NOTFOUND")
      find_path(BLAS_mkl.h_INCLUDE_DIRS
        NAMES mkl.h
        HINTS ${_inc_env})
    endif()
  endif()
  mark_as_advanced(BLAS_mkl.h_INCLUDE_DIRS)
  ## Print status if not found
  ## -------------------------
  #if (NOT BLAS_mkl.h_INCLUDE_DIRS AND MORSE_VERBOSE)
  #    Print_Find_Header_Status(blas mkl.h)
  #endif ()
  set(BLAS_INCLUDE_DIRS "")
  if(BLAS_mkl.h_INCLUDE_DIRS)
    list(APPEND BLAS_INCLUDE_DIRS "${BLAS_mkl.h_INCLUDE_DIRS}" )
  endif()

  ###
  # look for libs
  ###
  # if Intel 10 64 bit -> look for sequential and multithreaded versions
  if(BLA_VENDOR MATCHES "Intel10_64lp*")

    ## look for the sequential version
    set(BLA_VENDOR "Intel10_64lp_seq")
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "Look for the sequential version Intel10_64lp_seq")
    endif()
    find_package_blas(0)
    if(BLAS_FOUND)
      set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
      if (BLAS_CFLAGS_OTHER)
        set (BLAS_SEQ_CFLAGS_OTHER "${BLAS_CFLAGS_OTHER}")
      endif()
      if (BLAS_LDFLAGS_OTHER)
        set (BLAS_SEQ_LDFLAGS_OTHER "${BLAS_LDFLAGS_OTHER}")
      endif()
    else()
      set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
    endif()

    ## look for the multithreaded version
    set(BLA_VENDOR "Intel10_64lp")
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "Look for the multithreaded version Intel10_64lp")
    endif()
    find_package_blas(0)
    if(BLAS_FOUND)
      set(BLAS_PAR_LIBRARIES "${BLAS_LIBRARIES}")
      if (BLAS_CFLAGS_OTHER)
        set (BLAS_PAR_CFLAGS_OTHER "${BLAS_CFLAGS_OTHER}")
      endif()
      if (BLAS_LDFLAGS_OTHER)
        set (BLAS_PAR_LDFLAGS_OTHER "${BLAS_LDFLAGS_OTHER}")
      endif()
    else()
      set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")
    endif()

  else()

    if(BLAS_FOUND)
      set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
    else()
      set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
    endif()

  endif()

  # ACML case
elseif(BLA_VENDOR MATCHES "ACML*")

  ## look for the sequential version
  set(BLA_VENDOR "ACML")
  find_package_blas(0)
  if(BLAS_FOUND)
    set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
  else()
    set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
  endif()

  ## look for the multithreaded version
  set(BLA_VENDOR "ACML_MP")
  find_package_blas(0)
  if(BLAS_FOUND)
    set(BLAS_PAR_LIBRARIES "${BLAS_LIBRARIES}")
  else()
    set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")
  endif()

  # IBMESSL case
elseif(BLA_VENDOR MATCHES "IBMESSL*")

  ## look for the sequential version
  set(BLA_VENDOR "IBMESSL")
  find_package_blas(0)
  if(BLAS_FOUND)
    set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
  else()
    set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
  endif()

  ## look for the multithreaded version
  set(BLA_VENDOR "IBMESSLMT")
  find_package_blas(0)
  if(BLAS_FOUND)
    set(BLAS_PAR_LIBRARIES "${BLAS_LIBRARIES}")
  else()
    set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")
  endif()

else()

  if(BLAS_FOUND)
    # define the SEQ libs as the BLAS_LIBRARIES
    set(BLAS_SEQ_LIBRARIES "${BLAS_LIBRARIES}")
  else()
    set(BLAS_SEQ_LIBRARIES "${BLAS_SEQ_LIBRARIES-NOTFOUND}")
  endif()
  set(BLAS_PAR_LIBRARIES "${BLAS_PAR_LIBRARIES-NOTFOUND}")

endif()

# Reset pure BLAS cmake variables to the sequential case (arbitrary default)
if(BLAS_SEQ_LIBRARIES)
  set(BLAS_LIBRARIES "${BLAS_SEQ_LIBRARIES}")
  set(BLAS_CFLAGS_OTHER "${BLAS_SEQ_CFLAGS_OTHER}")
  set(BLAS_LDFLAGS_OTHER "${BLAS_SEQ_LDFLAGS_OTHER}")
endif()

# extract libs paths if not given by find_package(BLAS)
if (NOT BLAS_LIBRARY_DIRS)
  set(BLAS_LIBRARY_DIRS "")
  string(REPLACE " " ";" BLAS_LIBRARIES "${BLAS_LIBRARIES}")
  foreach(blas_lib ${BLAS_LIBRARIES})
    if (EXISTS "${blas_lib}")
      get_filename_component(a_blas_lib_dir "${blas_lib}" PATH)
      list(APPEND BLAS_LIBRARY_DIRS "${a_blas_lib_dir}" )
    else()
      string(REPLACE "-L" "" blas_lib "${blas_lib}")
      if (EXISTS "${blas_lib}")
        list(APPEND BLAS_LIBRARY_DIRS "${blas_lib}" )
      else()
        get_filename_component(a_blas_lib_dir "${blas_lib}" PATH)
        if (EXISTS "${a_blas_lib_dir}")
          list(APPEND BLAS_LIBRARY_DIRS "${a_blas_lib_dir}" )
        endif()
      endif()
    endif()
  endforeach()
  if (BLAS_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES BLAS_LIBRARY_DIRS)
  endif ()
endif(NOT BLAS_LIBRARY_DIRS)

# check that BLASEXT has been found
# ---------------------------------
include(FindPackageHandleStandardArgs)
if(BLA_VENDOR MATCHES "Intel*")
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "BLAS found is Intel MKL")
    message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(BLASEXT DEFAULT_MSG
    BLAS_SEQ_LIBRARIES
    BLAS_LIBRARY_DIRS
    BLAS_INCLUDE_DIRS)
  if(BLA_VENDOR MATCHES "Intel10_64lp*" AND BLAS_PAR_LIBRARIES)
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "BLAS parallel libraries stored in BLAS_PAR_LIBRARIES")
    endif()
    find_package_handle_standard_args(BLASEXT DEFAULT_MSG
      BLAS_PAR_LIBRARIES)
  endif()
elseif(BLA_VENDOR MATCHES "ACML*")
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "BLAS found is ACML")
    message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(BLASEXT DEFAULT_MSG
    BLAS_SEQ_LIBRARIES
    BLAS_LIBRARY_DIRS)
  if(BLAS_PAR_LIBRARIES)
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "BLAS parallel libraries stored in BLAS_PAR_LIBRARIES")
    endif()
    find_package_handle_standard_args(BLASEXT DEFAULT_MSG
      BLAS_PAR_LIBRARIES)
  endif()
elseif(BLA_VENDOR MATCHES "IBMESSL*")
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "BLAS found is ESSL")
    message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(BLASEXT DEFAULT_MSG
    BLAS_SEQ_LIBRARIES
    BLAS_LIBRARY_DIRS)
  if(BLAS_PAR_LIBRARIES)
    if(NOT BLASEXT_FIND_QUIETLY)
      message(STATUS "BLAS parallel libraries stored in BLAS_PAR_LIBRARIES")
    endif()
    find_package_handle_standard_args(BLASEXT DEFAULT_MSG
      BLAS_PAR_LIBRARIES)
  endif()
else()
  if(NOT BLASEXT_FIND_QUIETLY)
    message(STATUS "BLAS sequential libraries stored in BLAS_SEQ_LIBRARIES")
  endif()
  find_package_handle_standard_args(BLASEXT DEFAULT_MSG
    BLAS_SEQ_LIBRARIES
    BLAS_LIBRARY_DIRS)
endif()

if (BLAS_LIBRARIES)
  set(BLASEXT_LIBRARIES ${BLAS_LIBRARIES})
endif()
if (BLAS_INCLUDE_DIRS)
  set(BLASEXT_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
endif()
if (BLAS_LIBRARY_DIRS)
  set(BLASEXT_LIBRARY_DIRS ${BLAS_LIBRARY_DIRS})
endif()
