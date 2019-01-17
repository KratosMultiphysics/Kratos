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
# - Find TMG include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(TMG
#               [REQUIRED]             # Fail with error if tmg is not found
#              )
#
# This module finds headers and tmg library.
# Results are reported in variables:
#  TMG_FOUND            - True if headers and requested libraries were found
#  TMG_CFLAGS_OTHER      - tmglib compiler flags without headers paths
#  TMG_LDFLAGS_OTHER     - tmglib linker flags without libraries
#  TMG_INCLUDE_DIRS      - tmglib include directories
#  TMG_LIBRARY_DIRS      - tmglib link directories
#  TMG_LIBRARIES         - tmglib libraries to be linked (absolute path)
#  TMG_CFLAGS_OTHER_DEP  - tmglib + dependencies compiler flags without headers paths
#  TMG_LDFLAGS_OTHER_DEP - tmglib + dependencies linker flags without libraries
#  TMG_INCLUDE_DIRS_DEP  - tmglib + dependencies include directories
#  TMG_LIBRARY_DIRS_DEP  - tmglib + dependencies link directories
#  TMG_LIBRARIES_DEP     - tmglib + dependencies libraries
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DTMG=path/to/tmg):
#  TMG_DIR              - Where to find the base directory of tmg
#  TMG_INCDIR           - Where to find the header files
#  TMG_LIBDIR           - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: TMG_DIR, TMG_INCDIR, TMG_LIBDIR

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

if (NOT TMG_FOUND)
  set(TMG_DIR "" CACHE PATH "Installation directory of TMG library")
  if (NOT TMG_FIND_QUIETLY)
    message(STATUS "A cache variable, namely TMG_DIR, has been set to specify the install directory of TMG")
  endif()
endif()


# used to test a TMG function after
get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFunctionExists)
else (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFortranFunctionExists)
endif (NOT _LANGUAGES_ MATCHES Fortran)

# TMG depends on LAPACK anyway, try to find it
if(TMG_FIND_REQUIRED)
  find_package(LAPACK REQUIRED)
else()
  find_package(LAPACK)
endif()

# TMG depends on LAPACK
if (LAPACK_FOUND)

  # check if a tmg function exists in the LAPACK lib
  set(CMAKE_REQUIRED_LIBRARIES "${LAPACK_LDFLAGS_OTHER_DEP};${LAPACK_LIBRARIES_DEP}")
  set(CMAKE_REQUIRED_INCLUDES "${LAPACK_INCLUDE_DIRS_DEP}")
  set(CMAKE_REQUIRED_FLAGS "${LAPACK_CFLAGS_OTHER_DEP}")
  include(CheckFunctionExists)
  include(CheckFortranFunctionExists)
  unset(TMG_WORKS CACHE)
  if (NOT _LANGUAGES_ MATCHES Fortran)
    check_function_exists(dlarnv TMG_WORKS)
  else (NOT _LANGUAGES_ MATCHES Fortran)
    check_fortran_function_exists(dlarnv TMG_WORKS)
  endif (NOT _LANGUAGES_ MATCHES Fortran)
  if (TMG_WORKS)
    unset(TMG_WORKS CACHE)
    if (NOT _LANGUAGES_ MATCHES Fortran)
      check_function_exists(dlagsy TMG_WORKS)
    else (NOT _LANGUAGES_ MATCHES Fortran)
      check_fortran_function_exists(dlagsy TMG_WORKS)
    endif (NOT _LANGUAGES_ MATCHES Fortran)
    mark_as_advanced(TMG_WORKS)
  endif()

  if(TMG_WORKS)

    if(NOT TMG_FIND_QUIETLY)
      message(STATUS "Looking for tmg: test with lapack succeeds")
    endif()
    # test succeeds: TMG is in LAPACK
    set(TMG_LIBRARIES "${LAPACK_LIBRARIES}")

  else()

    if(NOT TMG_FIND_QUIETLY)
      message(STATUS "Looking for tmg : test with lapack fails")
      message(STATUS "Looking for tmg : try to find it elsewhere")
    endif()
    # test fails: try to find TMG lib exterior to LAPACK

    # Looking for lib tmg
    # -------------------

    # Add system library paths to search lib
    # --------------------------------------
    unset(_lib_env)
    set(ENV_TMG_DIR "$ENV{TMG_DIR}")
    set(ENV_TMG_LIBDIR "$ENV{TMG_LIBDIR}")
    if(ENV_TMG_LIBDIR)
      list(APPEND _lib_env "${ENV_TMG_LIBDIR}")
    elseif(ENV_TMG_DIR)
      list(APPEND _lib_env "${ENV_TMG_DIR}")
      list(APPEND _lib_env "${ENV_TMG_DIR}/lib")
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

    # Try to find the tmg lib in the given paths
    # ----------------------------------------------

    # call cmake macro to find the lib path
    if(TMG_LIBDIR)
      set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
      find_library(TMG_tmg_LIBRARY
        NAMES tmglib tmg
        HINTS ${TMG_LIBDIR} )
    else()
      if(TMG_DIR)
        set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
        find_library(TMG_tmg_LIBRARY
          NAMES tmglib tmg
          HINTS ${TMG_DIR}
          PATH_SUFFIXES lib lib32 lib64 )
      else()
        set(TMG_tmg_LIBRARY "TMG_tmg_LIBRARY-NOTFOUND")
        find_library(TMG_tmg_LIBRARY
          NAMES tmglib tmg
          HINTS ${_lib_env} )
      endif()
    endif()
    mark_as_advanced(TMG_tmg_LIBRARY)

    # If found, add path to cmake variable
    # ------------------------------------
    if (TMG_tmg_LIBRARY)
      get_filename_component(tmg_lib_path ${TMG_tmg_LIBRARY} PATH)
      # set cmake variables (respects naming convention)
      set(TMG_LIBRARIES    "${TMG_tmg_LIBRARY}")
      set(TMG_LIBRARY_DIRS "${tmg_lib_path}")
    else ()
      set(TMG_LIBRARIES    "TMG_LIBRARIES-NOTFOUND")
      set(TMG_LIBRARY_DIRS "TMG_LIBRARY_DIRS-NOTFOUND")
      if(NOT TMG_FIND_QUIETLY)
        message(STATUS "Looking for tmg -- lib tmg not found")
      endif()
    endif ()

    if (TMG_LIBRARY_DIRS)
      list(REMOVE_DUPLICATES TMG_LIBRARY_DIRS)
    endif ()

  endif(TMG_WORKS)

  # check a function to validate the find
  if(TMG_LIBRARIES)

    set(REQUIRED_INCDIRS)
    set(REQUIRED_FLAGS)
    set(REQUIRED_LDFLAGS)
    set(REQUIRED_LIBDIRS)
    set(REQUIRED_LIBS)

    # TMG
    if (TMG_INCLUDE_DIRS)
      set(REQUIRED_INCDIRS "${TMG_INCLUDE_DIRS}")
    endif()
    if (TMG_CFLAGS_OTHER)
      list(APPEND REQUIRED_FLAGS "${TMG_CFLAGS_OTHER}")
    endif()
    if (TMG_LDFLAGS_OTHER)
      list(APPEND REQUIRED_LDFLAGS "${TMG_LDFLAGS_OTHER}")
    endif()
    if (TMG_LIBRARY_DIRS)
      set(REQUIRED_LIBDIRS "${TMG_LIBRARY_DIRS}")
    endif()
    set(REQUIRED_LIBS "${TMG_LIBRARIES}")
    # LAPACK
    if (LAPACK_INCLUDE_DIRS_DEP)
      list(APPEND REQUIRED_INCDIRS "${LAPACK_INCLUDE_DIRS_DEP}")
    endif()
    if (LAPACK_CFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_FLAGS "${LAPACK_CFLAGS_OTHER_DEP}")
    endif()
    if (LAPACK_LDFLAGS_OTHER_DEP)
      list(APPEND REQUIRED_LDFLAGS "${LAPACK_LDFLAGS_OTHER_DEP}")
    endif()
    if (LAPACK_LIBRARY_DIRS_DEP)
      list(APPEND REQUIRED_LIBDIRS "${LAPACK_LIBRARY_DIRS_DEP}")
    endif()
    list(APPEND REQUIRED_LIBS "${LAPACK_LIBRARIES_DEP}")

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
    unset(TMG_WORKS CACHE)
    include(CheckFunctionExists)
    include(CheckFortranFunctionExists)
    if (NOT _LANGUAGES_ MATCHES Fortran)
      check_function_exists(dlarnv TMG_WORKS)
    else (NOT _LANGUAGES_ MATCHES Fortran)
      check_fortran_function_exists(dlarnv TMG_WORKS)
    endif (NOT _LANGUAGES_ MATCHES Fortran)
    if (TMG_WORKS)
      unset(TMG_WORKS CACHE)
      if (NOT _LANGUAGES_ MATCHES Fortran)
        check_function_exists(dlagsy TMG_WORKS)
      else (NOT _LANGUAGES_ MATCHES Fortran)
        check_fortran_function_exists(dlagsy TMG_WORKS)
      endif (NOT _LANGUAGES_ MATCHES Fortran)
      mark_as_advanced(TMG_WORKS)
    endif()

    if(TMG_WORKS)
      # save link with dependencies
      set(TMG_LIBRARIES_DEP "${REQUIRED_LIBS}")
      set(TMG_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
      set(TMG_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
      set(TMG_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
      set(TMG_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
    else()
      if(NOT TMG_FIND_QUIETLY)
        message(STATUS "Looking for tmg: test of dlarnv and dlagsy with tmg and lapack libraries fails")
        message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
        message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
        message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
        message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      endif()
    endif()
    set(CMAKE_REQUIRED_INCLUDES)
    set(CMAKE_REQUIRED_FLAGS)
    set(CMAKE_REQUIRED_LIBRARIES)
  endif(TMG_LIBRARIES)

else()

  if(NOT TMG_FIND_QUIETLY)
    message(STATUS "TMG requires LAPACK but LAPACK has not been found."
      "Please look for LAPACK first.")
  endif()

endif()

if (TMG_LIBRARIES)
  list(GET TMG_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT TMG_LIBRARY_DIRS)
    set(TMG_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(TMG_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of TMG library" FORCE)
  else()
    set(TMG_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of TMG library" FORCE)
  endif()
endif()
mark_as_advanced(TMG_DIR)
mark_as_advanced(TMG_DIR_FOUND)

# check that TMG has been found
# -------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TMG DEFAULT_MSG
  TMG_LIBRARIES
  TMG_WORKS)
