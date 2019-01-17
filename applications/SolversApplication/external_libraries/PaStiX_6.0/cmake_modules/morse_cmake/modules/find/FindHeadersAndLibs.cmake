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
#  @file FindHeadersAndLibs.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

# Some macros to print status when search for headers and libs
include(PrintFindStatus)

function(FindHeader _libname _header_to_find)

  # save _libname upper and lower case
  string(TOUPPER ${_libname} LIBNAME)
  string(TOLOWER ${_libname} libname)

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
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
  list(REMOVE_DUPLICATES _inc_env)


  # Try to find the _header_to_find in the given paths
  # --------------------------------------------------
  # call cmake macro to find the header path
  if(${LIBNAME}_INCDIR)
    set(${LIBNAME}_${_header_to_find}_DIRS "${LIBNAME}_${_header_to_find}_DIRS-NOTFOUND")
    find_path(${LIBNAME}_${_header_to_find}_DIRS
      NAMES ${_header_to_find}
      HINTS ${${LIBNAME}_INCDIR})
  elseif(${LIBNAME}_DIR)
    set(${LIBNAME}_${_header_to_find}_DIRS "${LIBNAME}_${_header_to_find}_DIRS-NOTFOUND")
    find_path(${LIBNAME}_${_header_to_find}_DIRS
      NAMES ${_header_to_find}
      HINTS ${${LIBNAME}_DIR}
      PATH_SUFFIXES include)
  else()
    set(${LIBNAME}_${_header_to_find}_DIRS "${LIBNAME}_${_header_to_find}_DIRS-NOTFOUND")
    find_path(${LIBNAME}_${_header_to_find}_DIRS
      NAMES ${_header_to_find}
      HINTS ${_inc_env})
  endif()
  mark_as_advanced(${LIBNAME}_${_header_to_find}_DIRS)

  # Print status if not found
  # -------------------------
  if (NOT ${LIBNAME}_${_header_to_find}_DIRS)
    Print_Find_Header_Status(${libname} ${_header_to_find})
  endif ()

endfunction(FindHeader)


##
## @end file FindHeadersAndLibs.cmake
##
