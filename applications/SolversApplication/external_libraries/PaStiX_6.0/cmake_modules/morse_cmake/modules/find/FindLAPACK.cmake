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
# - Find LAPACK library
# This module finds an installed fortran library that implements the LAPACK
# linear-algebra interface (see http://www.netlib.org/lapack/).
#
# The approach follows that taken for the autoconf macro file, acx_lapack.m4
# (distributed at http://ac-archive.sourceforge.net/ac-archive/acx_lapack.html).
#
# This module sets the following variables:
#  LAPACK_FOUND - set to true if a library implementing the LAPACK interface
#    is found
#  LAPACK_CFLAGS_OTHER      - lapack compiler flags without headers paths
#  LAPACK_LDFLAGS_OTHER     - lapack linker flags without libraries
#  LAPACK_INCLUDE_DIRS      - lapack include directories
#  LAPACK_LIBRARY_DIRS      - lapack link directories
#  LAPACK_LIBRARIES         - lapack libraries to be linked (absolute path)
#  LAPACK_LIBRARIES_SEQ     - lapack libraries to be linked (absolute path), sequential version (Intel MKL)
#  LAPACK_LIBRARIES_PAR     - lapack libraries to be linked (absolute path), multithreaded version (Intel MKL)
#  LAPACK_CFLAGS_OTHER_DEP  - lapack + dependencies compiler flags without headers paths
#  LAPACK_LDFLAGS_OTHER_DEP - lapack + dependencies linker flags without libraries
#  LAPACK_INCLUDE_DIRS_DEP  - lapack + dependencies include directories
#  LAPACK_LIBRARY_DIRS_DEP  - lapack + dependencies link directories
#  LAPACK_LIBRARIES_DEP     - lapack + dependencies libraries
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities
#  LAPACK_VENDOR_FOUND stores the LAPACK vendor found
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
#
#  LAPACK_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = LAPACK
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
# options at configure (ex: cmake path/to/project -DLAPACK_DIR=path/to/lapack):
#  LAPACK_DIR            - Where to find the base directory of lapack
#  LAPACK_INCDIR         - Where to find the header files
#  LAPACK_LIBDIR         - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: LAPACK_DIR, LAPACK_INCDIR, LAPACK_LIBDIR
# For MKL case and if no paths are given as hints, we will try to use the MKLROOT
# environment variable
# Note that if BLAS_DIR is set, it will also look for lapack in it
### List of vendors (BLA_VENDOR) valid in this module
##  Intel(mkl), ACML, Apple, NAS, Generic

#=============================================================================
# Copyright 2007-2009 Kitware, Inc.
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# Common macros to use in finds
include(FindInit)

if (NOT LAPACK_FOUND)
  set(LAPACK_DIR "" CACHE PATH "Installation directory of LAPACK library")
  if (NOT LAPACK_FIND_QUIETLY)
    message(STATUS "A cache variable, namely LAPACK_DIR, has been set to specify the install directory of LAPACK")
  endif()
endif (NOT LAPACK_FOUND)

option(LAPACK_VERBOSE "Print some additional information during LAPACK
libraries detection" OFF)
mark_as_advanced(LAPACK_VERBOSE)
if (BLAS_VERBOSE)
  set(LAPACK_VERBOSE ON)
endif ()
set(_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

get_property(_LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES)
if (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFunctionExists)
else (NOT _LANGUAGES_ MATCHES Fortran)
  include(CheckFortranFunctionExists)
endif (NOT _LANGUAGES_ MATCHES Fortran)

set(LAPACK_FOUND FALSE)
set(LAPACK95_FOUND FALSE)

# TODO: move this stuff to separate module

macro(Check_Lapack_Libraries LIBRARIES _prefix _name _flags _list _blas _threads)
  # This macro checks for the existence of the combination of fortran libraries
  # given by _list.  If the combination is found, this macro checks (using the
  # Check_Fortran_Function_Exists macro) whether can link against that library
  # combination using the name of a routine given by _name using the linker
  # flags given by _flags.  If the combination of libraries is found and passes
  # the link test, LIBRARIES is set to the list of complete library paths that
  # have been found.  Otherwise, LIBRARIES is set to FALSE.

  # N.B. _prefix is the prefix applied to the names of all cached variables that
  # are generated internally and marked advanced by this macro.
  set(_libdir ${ARGN})
  set(_libraries_work TRUE)
  set(${LIBRARIES})
  set(_combined_name)
  set(ENV_MKLROOT "$ENV{MKLROOT}")
  set(ENV_BLAS_DIR "$ENV{BLAS_DIR}")
  set(ENV_BLAS_LIBDIR "$ENV{BLAS_LIBDIR}")
  set(ENV_LAPACK_DIR "$ENV{LAPACK_DIR}")
  set(ENV_LAPACK_LIBDIR "$ENV{LAPACK_LIBDIR}")
  if (NOT _libdir)
    if (BLAS_LIBDIR)
      list(APPEND _libdir "${BLAS_LIBDIR}")
    elseif (BLAS_DIR)
      list(APPEND _libdir "${BLAS_DIR}")
      list(APPEND _libdir "${BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        list(APPEND _libdir "${BLAS_DIR}/lib64")
        list(APPEND _libdir "${BLAS_DIR}/lib/intel64")
      else()
        list(APPEND _libdir "${BLAS_DIR}/lib32")
        list(APPEND _libdir "${BLAS_DIR}/lib/ia32")
      endif()
    elseif(ENV_BLAS_LIBDIR)
      list(APPEND _libdir "${ENV_BLAS_LIBDIR}")
    elseif(ENV_BLAS_DIR)
      list(APPEND _libdir "${ENV_BLAS_DIR}")
      list(APPEND _libdir "${ENV_BLAS_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        list(APPEND _libdir "${ENV_BLAS_DIR}/lib64")
        list(APPEND _libdir "${ENV_BLAS_DIR}/lib/intel64")
      else()
        list(APPEND _libdir "${ENV_BLAS_DIR}/lib32")
        list(APPEND _libdir "${ENV_BLAS_DIR}/lib/ia32")
      endif()
    endif()
    if (LAPACK_LIBDIR)
      list(APPEND _libdir "${LAPACK_LIBDIR}")
    elseif (LAPACK_DIR)
      list(APPEND _libdir "${LAPACK_DIR}")
      list(APPEND _libdir "${LAPACK_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        list(APPEND _libdir "${LAPACK_DIR}/lib64")
        list(APPEND _libdir "${LAPACK_DIR}/lib/intel64")
      else()
        list(APPEND _libdir "${LAPACK_DIR}/lib32")
        list(APPEND _libdir "${LAPACK_DIR}/lib/ia32")
      endif()
    elseif(ENV_LAPACK_LIBDIR)
      list(APPEND _libdir "${ENV_LAPACK_LIBDIR}")
    elseif(ENV_LAPACK_DIR)
      list(APPEND _libdir "${ENV_LAPACK_DIR}")
      list(APPEND _libdir "${ENV_LAPACK_DIR}/lib")
      if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
        list(APPEND _libdir "${ENV_LAPACK_DIR}/lib64")
        list(APPEND _libdir "${ENV_LAPACK_DIR}/lib/intel64")
      else()
        list(APPEND _libdir "${ENV_LAPACK_DIR}/lib32")
        list(APPEND _libdir "${ENV_LAPACK_DIR}/lib/ia32")
      endif()
    else()
      if (ENV_MKLROOT)
        list(APPEND _libdir "${ENV_MKLROOT}/lib")
        if("${CMAKE_SIZEOF_VOID_P}" EQUAL "8")
          list(APPEND _libdir "${ENV_MKLROOT}/lib64")
          list(APPEND _libdir "${ENV_MKLROOT}/lib/intel64")
        else()
          list(APPEND _libdir "${ENV_MKLROOT}/lib32")
          list(APPEND _libdir "${ENV_MKLROOT}/lib/ia32")
        endif()
      endif()
      if (WIN32)
        string(REPLACE ":" ";" _libdir2 "$ENV{LIB}")
      elseif (APPLE)
        string(REPLACE ":" ";" _libdir2 "$ENV{DYLD_LIBRARY_PATH}")
      else ()
        string(REPLACE ":" ";" _libdir2 "$ENV{LD_LIBRARY_PATH}")
      endif ()
      list(APPEND _libdir "${_libdir2}")
      list(APPEND _libdir "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
    endif ()
  endif ()

  if (LAPACK_VERBOSE)
    message("${Cyan}Try to find LAPACK libraries: ${_list}")
  endif ()

  foreach(_library ${_list})
    set(_combined_name ${_combined_name}_${_library})

    if(_libraries_work)
      if (BLA_STATIC)
        if (WIN32)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif ( WIN32 )
        if (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .lib ${CMAKE_FIND_LIBRARY_SUFFIXES})
        else (APPLE)
          set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
        endif (APPLE)
      else (BLA_STATIC)
        if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
          # for ubuntu's libblas3gf and liblapack3gf packages
          set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} .so.3gf)
        endif ()
      endif (BLA_STATIC)
      find_library(${_prefix}_${_library}_LIBRARY
        NAMES ${_library}
        HINTS ${_libdir}
        )
      mark_as_advanced(${_prefix}_${_library}_LIBRARY)
      # Print status if not found
      # -------------------------
      if (NOT ${_prefix}_${_library}_LIBRARY AND NOT LAPACK_FIND_QUIETLY AND LAPACK_VERBOSE)
        Print_Find_Library_Blas_Status(lapack ${_library} ${_libdir})
      endif ()
      set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
      set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
    endif(_libraries_work)
  endforeach(_library ${_list})

  if(_libraries_work)
    # Test this combination of libraries.
    if (CMAKE_SYSTEM_NAME STREQUAL "Linux" AND BLA_STATIC)
      list(INSERT ${LIBRARIES} 0 "-Wl,--start-group")
      list(APPEND ${LIBRARIES} "-Wl,--end-group")
    endif()
    if(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} "-Wl,--start-group" ${${LIBRARIES}} ${_blas} "-Wl,--end-group" ${_threads})
    else(UNIX AND BLA_STATIC)
      set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_blas} ${_threads})
    endif(UNIX AND BLA_STATIC)
    set(CMAKE_REQUIRED_INCLUDES "${LAPACK_INCLUDE_DIRS}")
    set(CMAKE_REQUIRED_FLAGS "${LAPACK_CFLAGS_OTHER}")
    if (LAPACK_VERBOSE)
      message("${Cyan}LAPACK libs found. Try to compile symbol ${_name} with"
        "following libraries: ${CMAKE_REQUIRED_LIBRARIES}")
    endif ()
    if(NOT LAPACK_FOUND)
      unset(${_prefix}${_combined_name}_WORKS CACHE)
    endif()
    if (NOT _LANGUAGES_ MATCHES Fortran)
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    else (NOT _LANGUAGES_ MATCHES Fortran)
      check_fortran_function_exists(${_name} ${_prefix}${_combined_name}_WORKS)
    endif (NOT _LANGUAGES_ MATCHES Fortran)
    set(CMAKE_REQUIRED_LIBRARIES)
    mark_as_advanced(${_prefix}${_combined_name}_WORKS)
    set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
  endif(_libraries_work)

  if(_libraries_work)
    set(${LIBRARIES} ${${LIBRARIES}} ${_blas} ${_threads})
  else(_libraries_work)
    set(${LIBRARIES} FALSE)
  endif(_libraries_work)

endmacro(Check_Lapack_Libraries)


set(LAPACK_LDFLAGS_OTHER)
set(LAPACK_LIBRARIES)
set(LAPACK95_LIBRARIES)

if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
  find_package(BLAS)
else(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
  find_package(BLAS REQUIRED)
endif(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)

if(BLAS_FOUND)
  set(LAPACK_INCLUDE_DIRS ${BLAS_INCLUDE_DIRS})
  if ($ENV{BLA_VENDOR} MATCHES ".+")
    set(BLA_VENDOR $ENV{BLA_VENDOR})
  else ($ENV{BLA_VENDOR} MATCHES ".+")
    if(NOT BLA_VENDOR)
      set(BLA_VENDOR "All")
    endif(NOT BLA_VENDOR)
  endif ($ENV{BLA_VENDOR} MATCHES ".+")

  if (UNIX AND NOT WIN32)
    # m
    find_library(M_LIBRARY NAMES m)
    mark_as_advanced(M_LIBRARY)
    if(M_LIBRARY)
      set(LM "${M_LIBRARY}")
    else()
      set(LM "")
    endif()
  endif()

  # check if blas lib contains lapack symbols
  if (BLAS_LIBRARIES)
    if (BLAS_LIBRARIES_DEP)
      set(LIBRARIES ${BLAS_LIBRARIES_DEP})
    else()
      set(LIBRARIES ${BLAS_LIBRARIES})
    endif()
    check_lapack_libraries(
      LAPACK_LIBRARIES
      LAPACK
      "CHEEV"
      ""
      ""
      "${LIBRARIES}"
      ""
      )
    if (LAPACK_LIBRARIES)
      if (BLAS_LIBRARIES_SEQ)
        set(LAPACK_LIBRARIES_SEQ "${BLAS_LIBRARIES_SEQ}")
      endif()
      if (BLAS_LIBRARIES_PAR)
        set(LAPACK_LIBRARIES_PAR "${BLAS_LIBRARIES_PAR}")
      endif()
      if(NOT LAPACK_FIND_QUIETLY)
            if(LAPACK_LIBRARIES)
              message(STATUS "Looking for LAPACK in BLAS: found")
            else()
              message(STATUS "Looking for LAPACK in BLAS: not found")
            endif()
          endif()
      if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
          set (LAPACK_VENDOR_FOUND "${BLAS_VENDOR_FOUND}")
      endif()
    endif()
  endif(BLAS_LIBRARIES)

  # if not lapack in blas libs, try to find lapack with pkg-config
  set(ENV_LAPACK_DIR "$ENV{LAPACK_DIR}")
  set(ENV_MKL_DIR "$ENV{MKLROOT}")
  set(ENV_LAPACK_INCDIR "$ENV{LAPACK_INCDIR}")
  set(ENV_LAPACK_LIBDIR "$ENV{LAPACK_LIBDIR}")
  set(LAPACK_GIVEN_BY_USER "FALSE")
  if ( LAPACK_DIR OR ( LAPACK_INCDIR AND LAPACK_LIBDIR) OR ENV_LAPACK_DIR OR ENV_MKL_DIR OR (ENV_LAPACK_INCDIR AND ENV_LAPACK_LIBDIR) )
    set(LAPACK_GIVEN_BY_USER "TRUE")
  endif()

  # Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
  # -------------------------------------------------------------------------------------
  include(FindPkgConfig)
  find_package(PkgConfig QUIET)
  if( PKG_CONFIG_EXECUTABLE AND NOT LAPACK_GIVEN_BY_USER AND BLA_VENDOR STREQUAL "All")

    if (BLA_STATIC)
      set(MKL_STR_BLA_STATIC "static")
    else()
      set(MKL_STR_BLA_STATIC "dynamic")
    endif()
    # try different blas
    if (BLA_VENDOR STREQUAL "Intel10_64lp")
      pkg_search_module(LAPACK mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Intel MKL Sequential")
        set (LAPACK_LIBRARIES_SEQ "${LAPACK_LIBRARIES}")
      endif()
      pkg_search_module(LAPACK mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Intel MKL Multithreaded")
        set (LAPACK_LIBRARIES_PAR "${LAPACK_LIBRARIES}")
      endif()
    elseif(BLA_VENDOR STREQUAL "Intel10_64lp_seq")
      pkg_search_module(LAPACK mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Intel MKL")
        set (LAPACK_LIBRARIES_SEQ "${LAPACK_LIBRARIES}")
      endif()
    elseif(BLA_VENDOR STREQUAL "Open")
      pkg_search_module(LAPACK openblas)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Openblas")
      endif()
    elseif(BLA_VENDOR STREQUAL "Generic")
      pkg_search_module(LAPACK lapack)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Netlib or other Generic liblapack")
      endif()
    else()
      if (NOT LAPACK_FOUND)
      pkg_search_module(LAPACK mkl-${MKL_STR_BLA_STATIC}-lp64-iomp)
      if (LAPACK_FOUND)
        set (LAPACK_VENDOR_FOUND "Intel MKL Multithreaded")
        set (LAPACK_LIBRARIES_PAR "${LAPACK_LIBRARIES}")
      endif()
        pkg_search_module(LAPACK mkl-${MKL_STR_BLA_STATIC}-lp64-seq)
        if (LAPACK_FOUND)
          set (LAPACK_VENDOR_FOUND "Intel MKL")
          set (LAPACK_LIBRARIES_SEQ "${LAPACK_LIBRARIES}")
        endif()
      endif()
      if (NOT LAPACK_FOUND)
        pkg_search_module(LAPACK openblas)
        if (LAPACK_FOUND)
          set (LAPACK_VENDOR_FOUND "Openblas")
        endif()
      endif()
      if (NOT LAPACK_FOUND)
        pkg_search_module(LAPACK lapack)
        if (LAPACK_FOUND)
          set (LAPACK_VENDOR_FOUND "Netlib or other Generic liblapack")
        endif()
      endif()
    endif()

    if (NOT LAPACK_FIND_QUIETLY)
      if (LAPACK_FOUND AND LAPACK_LIBRARIES)
        message(STATUS "Looking for LAPACK - found using PkgConfig")
      else()
        message(STATUS "${Magenta}Looking for LAPACK - not found using PkgConfig."
          "\n   Perhaps you should add the directory containing lapack.pc to"
          "\n   the PKG_CONFIG_PATH environment variable.${ColourReset}")
      endif()
    endif()

    if (LAPACK_FOUND AND LAPACK_LIBRARIES)
      set(LAPACK_FOUND_WITH_PKGCONFIG "TRUE")
      find_pkgconfig_libraries_absolute_path(LAPACK)
    else()
      set(LAPACK_FOUND_WITH_PKGCONFIG "FALSE")
    endif()

  endif()

  if( (NOT LAPACK_FOUND_WITH_PKGCONFIG) OR LAPACK_GIVEN_BY_USER )
    #intel lapack
    if (BLA_VENDOR MATCHES "Intel" OR BLA_VENDOR STREQUAL "All")

      if(NOT LAPACK_LIBRARIES)
        if (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)

          if(LAPACK_FIND_QUIETLY OR NOT LAPACK_FIND_REQUIRED)
                find_PACKAGE(Threads)
          else()
                find_package(Threads REQUIRED)
          endif()
          if( THREADS_FOUND )
            libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
          endif ()

          set(LAPACK_SEARCH_LIBS "")

          set(additional_flags "")
          if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
                set(additional_flags "-Wl,--no-as-needed")
          endif()

          if (BLA_F95)
                set(LAPACK_mkl_SEARCH_SYMBOL "CHEEV")
                set(_LIBRARIES LAPACK95_LIBRARIES)
                set(_BLAS_LIBRARIES ${BLAS95_LIBRARIES})

                # old
                list(APPEND LAPACK_SEARCH_LIBS
                  "mkl_lapack95")
                # new >= 10.3
                list(APPEND LAPACK_SEARCH_LIBS
                  "mkl_intel_c")
                list(APPEND LAPACK_SEARCH_LIBS
                  "mkl_intel_lp64")
          else(BLA_F95)
                set(LAPACK_mkl_SEARCH_SYMBOL "cheev")
                set(_LIBRARIES LAPACK_LIBRARIES)
                set(_BLAS_LIBRARIES ${BLAS_LIBRARIES})

                # old
                list(APPEND LAPACK_SEARCH_LIBS
                  "mkl_lapack")
                # new >= 10.3
                list(APPEND LAPACK_SEARCH_LIBS
                  "mkl_gf_lp64")
          endif(BLA_F95)

          # First try empty lapack libs
          if (NOT ${_LIBRARIES})
                check_lapack_libraries(
                  ${_LIBRARIES}
                  LAPACK
                  ${LAPACK_mkl_SEARCH_SYMBOL}
                  "${additional_flags}"
                  ""
                  "${_BLAS_LIBRARIES}"
                  "${CMAKE_THREAD_LIBS_INIT};${LM}"
                  )
                if(_LIBRARIES)
                  set(LAPACK_LDFLAGS_OTHER "${additional_flags}")
                endif()
          endif ()
          # Then try the search libs
          foreach (IT ${LAPACK_SEARCH_LIBS})
                if (NOT ${_LIBRARIES})
                  check_lapack_libraries(
                    ${_LIBRARIES}
                    LAPACK
                    ${LAPACK_mkl_SEARCH_SYMBOL}
                    "${additional_flags}"
                    "${IT}"
                    "${_BLAS_LIBRARIES}"
                    "${CMAKE_THREAD_LIBS_INIT};${LM}"
                    )
                  if(_LIBRARIES)
                    set(LAPACK_LDFLAGS_OTHER "${additional_flags}")
                  endif()
                endif ()
          endforeach ()
          if(NOT LAPACK_FIND_QUIETLY)
                if(${_LIBRARIES})
                  message(STATUS "Looking for MKL LAPACK: found")
                else()
                  message(STATUS "Looking for MKL LAPACK: not found")
                endif()
          endif(NOT LAPACK_FIND_QUIETLY)
          if(${_LIBRARIES} AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "Intel MKL")
          endif()

        endif (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
      endif(NOT LAPACK_LIBRARIES)
    endif(BLA_VENDOR MATCHES "Intel" OR BLA_VENDOR STREQUAL "All")

    #goto lapack
    if (BLA_VENDOR STREQUAL "Goto" OR BLA_VENDOR STREQUAL "All")
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheev
          ""
          "goto2"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for Goto LAPACK: found")
          else()
            message(STATUS "Looking for Goto LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "Goto")
        endif()
      endif(NOT LAPACK_LIBRARIES)
    endif (BLA_VENDOR STREQUAL "Goto" OR BLA_VENDOR STREQUAL "All")

    #open lapack
    if (BLA_VENDOR STREQUAL "Open" OR BLA_VENDOR STREQUAL "All")
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheev
          ""
          "openblas"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for Open LAPACK: found")
          else()
            message(STATUS "Looking for Open LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "Openblas")
        endif()
      endif(NOT LAPACK_LIBRARIES)
    endif (BLA_VENDOR STREQUAL "Open" OR BLA_VENDOR STREQUAL "All")

    # LAPACK in IBM ESSL library (requires generic lapack lib, too)
    if (BLA_VENDOR STREQUAL "IBMESSL" OR BLA_VENDOR STREQUAL "All")
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheevd
          ""
          "essl"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for IBM ESSL LAPACK: found")
          else()
            message(STATUS "Looking for IBM ESSL LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "IBM ESSL")
        endif()
      endif()
    endif ()

    # LAPACK in IBM ESSL_MT library (requires generic lapack lib, too)
    if (BLA_VENDOR STREQUAL "IBMESSLMT" OR BLA_VENDOR STREQUAL "All")
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheevd
          ""
          "esslsmp"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for IBM ESSL MT LAPACK: found")
          else()
            message(STATUS "Looking for IBM ESSL MT LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "IBM ESSL MT")
        endif()
      endif()
    endif ()

    #acml lapack
    if (BLA_VENDOR MATCHES "ACML.*" OR BLA_VENDOR STREQUAL "All")
      if (BLAS_LIBRARIES MATCHES ".+acml.+")
        set (LAPACK_LIBRARIES ${BLAS_LIBRARIES})
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for ACML LAPACK: found")
          else()
            message(STATUS "Looking for ACML LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "ACML")
        endif()
      endif ()
    endif ()

    # Apple LAPACK library?
    if (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")
      if(NOT LAPACK_LIBRARIES)
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheev
          ""
          "Accelerate"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for Apple Accelerate LAPACK: found")
          else()
            message(STATUS "Looking for Apple Accelerate LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "Apple Accelerate")
        endif()
      endif(NOT LAPACK_LIBRARIES)
    endif (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")

    if (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")
      if ( NOT LAPACK_LIBRARIES )
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheev
          ""
          "vecLib"
          "${BLAS_LIBRARIES}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for NAS LAPACK: found")
          else()
            message(STATUS "Looking for NAS LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "NAS")
        endif()
      endif ( NOT LAPACK_LIBRARIES )
    endif (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")

    # Generic LAPACK library?
    if (BLA_VENDOR STREQUAL "Generic" OR
        BLA_VENDOR STREQUAL "ATLAS" OR
        BLA_VENDOR STREQUAL "All")
      if ( NOT LAPACK_LIBRARIES )
        check_lapack_libraries(
          LAPACK_LIBRARIES
          LAPACK
          cheev
          ""
          "lapack"
          "${BLAS_LIBRARIES};${LM}"
          ""
          )
        if(NOT LAPACK_FIND_QUIETLY)
          if(LAPACK_LIBRARIES)
            message(STATUS "Looking for Generic LAPACK: found")
          else()
            message(STATUS "Looking for Generic LAPACK: not found")
          endif()
        endif()
        if (LAPACK_LIBRARIES AND NOT LAPACK_VENDOR_FOUND)
            set (LAPACK_VENDOR_FOUND "Netlib or other Generic liblapack")
        endif()
      endif ( NOT LAPACK_LIBRARIES )
    endif ()
  endif( (NOT LAPACK_FOUND_WITH_PKGCONFIG) OR LAPACK_GIVEN_BY_USER )

else(BLAS_FOUND)
  message(STATUS "LAPACK requires BLAS")
endif(BLAS_FOUND)

if(BLA_F95)
  if(LAPACK95_LIBRARIES)
    set(LAPACK95_FOUND TRUE)
  else(LAPACK95_LIBRARIES)
    set(LAPACK95_FOUND FALSE)
  endif(LAPACK95_LIBRARIES)
  if(NOT LAPACK_FIND_QUIETLY)
    if(LAPACK95_FOUND)
      message(STATUS "A library with LAPACK95 API found.")
      message(STATUS "LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
    else(LAPACK95_FOUND)
      message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but LAPACK 95 libraries could not be found or check of symbols failed."
        "\nPlease indicate where to find LAPACK libraries. You have three options:\n"
        "- Option 1: Provide the installation directory of LAPACK library with cmake option: -DLAPACK_DIR=your/path/to/lapack\n"
        "- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
        "- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
        "\nTo follow libraries detection more precisely you can activate a verbose mode with -DLAPACK_VERBOSE=ON at cmake configure."
        "\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
        "\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
        "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
        "Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
      if(LAPACK_FIND_REQUIRED)
        message(FATAL_ERROR
          "A required library with LAPACK95 API not found. Please specify library location."
          )
      else(LAPACK_FIND_REQUIRED)
        message(STATUS
          "A library with LAPACK95 API not found. Please specify library location."
          )
      endif(LAPACK_FIND_REQUIRED)
    endif(LAPACK95_FOUND)
  endif(NOT LAPACK_FIND_QUIETLY)
  set(LAPACK_FOUND "${LAPACK95_FOUND}")
  set(LAPACK_LIBRARIES "${LAPACK95_LIBRARIES}")
else(BLA_F95)
  if(LAPACK_LIBRARIES)
    set(LAPACK_FOUND TRUE)
  else(LAPACK_LIBRARIES)
    set(LAPACK_FOUND FALSE)
  endif(LAPACK_LIBRARIES)

  if(NOT LAPACK_FIND_QUIETLY)
    if(LAPACK_FOUND)
      message(STATUS "A library with LAPACK API found.")
      message(STATUS "LAPACK_LIBRARIES ${LAPACK_LIBRARIES}")
    else(LAPACK_FOUND)
      message(WARNING "BLA_VENDOR has been set to ${BLA_VENDOR} but LAPACK libraries could not be found or check of symbols failed."
        "\nPlease indicate where to find LAPACK libraries. You have three options:\n"
        "- Option 1: Provide the installation directory of LAPACK library with cmake option: -DLAPACK_DIR=your/path/to/lapack\n"
        "- Option 2: Provide the directory where to find BLAS libraries with cmake option: -DBLAS_LIBDIR=your/path/to/blas/libs\n"
        "- Option 3: Update your environment variable (Linux: LD_LIBRARY_PATH, Windows: LIB, Mac: DYLD_LIBRARY_PATH)\n"
        "\nTo follow libraries detection more precisely you can activate a verbose mode with -DLAPACK_VERBOSE=ON at cmake configure."
        "\nYou could also specify a BLAS vendor to look for by setting -DBLA_VENDOR=blas_vendor_name."
        "\nList of possible BLAS vendor: Goto, ATLAS PhiPACK, CXML, DXML, SunPerf, SCSL, SGIMATH, IBMESSL, Intel10_32 (intel mkl v10 32 bit),"
        "Intel10_64lp (intel mkl v10 64 bit, lp thread model, lp64 model), Intel10_64lp_seq (intel mkl v10 64 bit, sequential code, lp64 model),"
        "Intel( older versions of mkl 32 and 64 bit), ACML, ACML_MP, ACML_GPU, Apple, NAS, Generic")
      if(LAPACK_FIND_REQUIRED)
        message(FATAL_ERROR
          "A required library with LAPACK API not found. Please specify library location."
          )
      else(LAPACK_FIND_REQUIRED)
        message(STATUS
          "A library with LAPACK API not found. Please specify library location."
          )
      endif(LAPACK_FIND_REQUIRED)
    endif(LAPACK_FOUND)
  endif(NOT LAPACK_FIND_QUIETLY)
endif(BLA_F95)

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_lapack_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

if (LAPACK_FOUND)
  list(GET LAPACK_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (NOT LAPACK_LIBRARY_DIRS)
    set(LAPACK_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)")
    string(REGEX REPLACE "(/lib(32|64)?$)|(/lib/intel64$|/lib/ia32$)" "" not_cached_dir "${first_lib_path}")
    set(LAPACK_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of LAPACK library" FORCE)
  else()
    set(LAPACK_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of LAPACK library" FORCE)
  endif()
  set(LAPACK_CFLAGS_OTHER_DEP "${LAPACK_CFLAGS_OTHER}" "${BLAS_CFLAGS_OTHER}")
  set(LAPACK_LDFLAGS_OTHER_DEP "${LAPACK_LDFLAGS_OTHER}" "${BLAS_LDFLAGS_OTHER}")
  set(LAPACK_INCLUDE_DIRS_DEP "${LAPACK_INCLUDE_DIRS}" "${BLAS_INCLUDE_DIRS}")
  set(LAPACK_LIBRARY_DIRS_DEP "${LAPACK_LIBRARY_DIRS}" "${BLAS_LIBRARY_DIRS}")
  set(LAPACK_LIBRARIES_DEP "${LAPACK_LIBRARIES}" "${BLAS_LIBRARIES}")
  list(REMOVE_DUPLICATES LAPACK_CFLAGS_OTHER_DEP)
  list(REMOVE_DUPLICATES LAPACK_LDFLAGS_OTHER_DEP)
  list(REMOVE_DUPLICATES LAPACK_INCLUDE_DIRS_DEP)
  list(REMOVE_DUPLICATES LAPACK_LIBRARY_DIRS_DEP)
endif()
mark_as_advanced(LAPACK_DIR)
mark_as_advanced(LAPACK_DIR_FOUND)
