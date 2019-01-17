###
#
# @copyright (c) 2009-2015 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2018 Inria. All rights reserved.
# @copyright (c) 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
# - Find PARSEC include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(PARSEC
#               [version] [EXACT]      # Minimum or EXACT version e.g. 1.1
#               [REQUIRED]             # Fail with error if parsec is not found
#               [COMPONENTS <comp1> <comp2> ...] # dependencies
#              )
#
#  PARSEC depends on the following libraries:
#   - Threads, m, rt
#
#  COMPONENTS are optional libraries PARSEC could be linked with,
#  Use it to drive detection of a specific compilation chain
#  COMPONENTS can be some of the following:
#   - HWLOC: to activate the detection of PARSEC linked with HWLOC
#   - CUDA: to activate the detection of PARSEC linked with CUDA
#   - MPI: to activate the detection of PARSEC linked with MPI
#   - AYUDAME: ??
#
# Results are reported in variables:
#  PARSEC_FOUND             - True if headers and requested libraries were found
#  PARSEC_CFLAGS_OTHER      - parsec compiler flags without headers paths
#  PARSEC_LDFLAGS_OTHER     - parsec linker flags without libraries
#  PARSEC_INCLUDE_DIRS      - parsec include directories
#  PARSEC_LIBRARY_DIRS      - parsec link directories
#  PARSEC_LIBRARIES         - parsec libraries to be linked (absolute path)
#  PARSEC_CFLAGS_OTHER_DEP  - parsec + dependencies compiler flags without headers paths
#  PARSEC_LDFLAGS_OTHER_DEP - parsec + dependencies linker flags without libraries
#  PARSEC_INCLUDE_DIRS_DEP  - parsec + dependencies include directories
#  PARSEC_LIBRARY_DIRS_DEP  - parsec + dependencies link directories
#  PARSEC_LIBRARIES_DEP     - parsec + dependencies libraries
#
#  PARSEC_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#  if found with pkg-config the following variables are set
#  <PREFIX>  = PARSEC
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
# options at configure (ex: cmake path/to/project -DPARSEC=path/to/parsec):
#  PARSEC_DIR                    - Where to find the base directory of parsec
#  PARSEC_INCDIR                 - Where to find the header files
#  PARSEC_LIBDIR                 - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: PARSEC_DIR, PARSEC_INCDIR, PARSEC_LIBDIR

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

include(CheckSymbolExists)

if (NOT PARSEC_FOUND)
  set(PARSEC_DIR "" CACHE PATH "Installation directory of PARSEC library")
  if (NOT PARSEC_FIND_QUIETLY)
    message(STATUS "A cache variable, namely PARSEC_DIR, has been set to specify the install directory of PARSEC")
  endif()
endif()

# PARSEC may depend on other packages (HWLOC, MPI, CUDA, ...)
# try to find them if specified as COMPONENTS during the call
set(PARSEC_LOOK_FOR_HWLOC FALSE)
set(PARSEC_LOOK_FOR_CUDA FALSE)
set(PARSEC_LOOK_FOR_MPI FALSE)

if( PARSEC_FIND_COMPONENTS )
  foreach( component ${PARSEC_FIND_COMPONENTS} )
    if(${component} STREQUAL "HWLOC")
      set(PARSEC_LOOK_FOR_HWLOC TRUE)
    elseif(${component} STREQUAL "CUDA")
      set(PARSEC_LOOK_FOR_CUDA TRUE)
    elseif(${component} STREQUAL "MPI")
      set(PARSEC_LOOK_FOR_MPI TRUE)
    endif()
  endforeach()
endif()

# Required dependencies
# ---------------------

if (NOT PARSEC_FIND_QUIETLY)
  message(STATUS "Looking for PARSEC - Try to detect pthread")
endif()
if (PARSEC_FIND_REQUIRED)
  find_package(Threads REQUIRED)
else()
  find_package(Threads)
endif()
if( THREADS_FOUND )
  libraries_absolute_path(CMAKE_THREAD_LIBS_INIT "")
endif ()
set(PARSEC_EXTRA_LIBRARIES "")
if( THREADS_FOUND )
  list(APPEND PARSEC_EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
if (NOT PARSEC_FIND_QUIETLY)
  message(STATUS "Looking for PARSEC - Try to detect libm")
endif()
set(PARSEC_M_LIBRARIES "")
if(UNIX OR WIN32)
  find_library(
    PARSEC_M_m_LIBRARY
    NAMES m
    )
  mark_as_advanced(PARSEC_M_m_LIBRARY)
  if (PARSEC_M_m_LIBRARY)
    list(APPEND PARSEC_M_LIBRARIES "${PARSEC_M_m_LIBRARY}")
    list(APPEND PARSEC_EXTRA_LIBRARIES "${PARSEC_M_m_LIBRARY}")
  else()
    if (PARSEC_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find libm on your system."
        "Are you sure to a have a C compiler installed?")
    endif()
  endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
if (NOT PARSEC_FIND_QUIETLY)
  message(STATUS "Looking for PARSEC - Try to detect librt")
endif()
set(PARSEC_RT_LIBRARIES "")
if(UNIX AND NOT APPLE)
  find_library(
    PARSEC_RT_rt_LIBRARY
    NAMES rt
    )
  mark_as_advanced(PARSEC_RT_rt_LIBRARY)
  if (PARSEC_RT_rt_LIBRARY)
    list(APPEND PARSEC_RT_LIBRARIES "${PARSEC_RT_rt_LIBRARY}")
    list(APPEND PARSEC_EXTRA_LIBRARIES "${PARSEC_RT_rt_LIBRARY}")
  else()
    if (PARSEC_FIND_REQUIRED)
      message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
  endif()
endif()

# Try to find libdl
if (NOT PARSEC_FIND_QUIETLY)
  message(STATUS "Looking for PARSEC - Try to detect libdl")
endif()
find_library(
  DL_LIBRARY
  NAMES dl
  )
mark_as_advanced(DL_LIBRARY)
if (DL_LIBRARY)
  list(APPEND PARSEC_EXTRA_LIBRARIES "${DL_LIBRARY}")
endif()

# PARSEC may depend on HWLOC, try to find it
if (PARSEC_LOOK_FOR_HWLOC)
  if (PARSEC_FIND_REQUIRED)
    find_package(HWLOC REQUIRED)
  else()
    find_package(HWLOC)
  endif()
endif()

# PARSEC may depend on CUDA, try to find it
if (PARSEC_LOOK_FOR_CUDA)
  if (PARSEC_FIND_REQUIRED AND PARSEC_FIND_REQUIRED_CUDA)
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

# PARSEC may depend on MPI, try to find it
if (PARSEC_LOOK_FOR_MPI)
  if (PARSEC_FIND_REQUIRED AND PARSEC_FIND_REQUIRED_MPI)
    find_package(MPI REQUIRED)
  else()
    find_package(MPI)
  endif()
  if (MPI_FOUND)
    mark_as_advanced(MPI_LIBRARY)
    mark_as_advanced(MPI_EXTRA_LIBRARY)
  endif()
endif()

set(ENV_PARSEC_DIR "$ENV{PARSEC_DIR}")
set(ENV_PARSEC_INCDIR "$ENV{PARSEC_INCDIR}")
set(ENV_PARSEC_LIBDIR "$ENV{PARSEC_LIBDIR}")
set(PARSEC_GIVEN_BY_USER "FALSE")
if ( PARSEC_DIR OR ( PARSEC_INCDIR AND PARSEC_LIBDIR) OR ENV_PARSEC_DIR OR (ENV_PARSEC_INCDIR AND ENV_PARSEC_LIBDIR) )
  set(PARSEC_GIVEN_BY_USER "TRUE")
endif()

# Optionally use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)

if(PKG_CONFIG_EXECUTABLE AND NOT PARSEC_GIVEN_BY_USER)

  pkg_search_module(PARSEC parsec)

  if (NOT PARSEC_FIND_QUIETLY)
    if (PARSEC_FOUND AND PARSEC_LIBRARIES)
      message(STATUS "Looking for PARSEC - found using PkgConfig")
    else()
      message(STATUS "${Magenta}Looking for PARSEC - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing libparsec.pc"
        "\n   to the PKG_CONFIG_PATH environment variable.${ColourReset}")
    endif()
  endif()

  # create list of binaries to find
  set(PARSEC_bins_to_find "parsec_ptgpp")

  # call cmake macro to find the bin path
  if(PARSEC_PREFIX)
    foreach(parsec_bin ${PARSEC_bins_to_find})
      set(PARSEC_${parsec_bin}_BIN_DIR "PARSEC_${parsec_bin}_BIN_DIR-NOTFOUND")
      find_path(PARSEC_${parsec_bin}_BIN_DIR
        NAMES ${parsec_bin}
        HINTS ${PARSEC_PREFIX}
        PATH_SUFFIXES "bin")
    endforeach()
  else()
    if (PARSEC_FIND_REQUIRED)
      message(FATAL_ERROR "PARSEC_PREFIX not defined by pkg_search_module")
    endif()
  endif()

  if (PARSEC_FOUND AND PARSEC_LIBRARIES)
    set(PARSEC_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(PARSEC)
  else()
    set(PARSEC_FOUND_WITH_PKGCONFIG "FALSE")
  endif()

endif(PKG_CONFIG_EXECUTABLE AND NOT PARSEC_GIVEN_BY_USER)

if( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PARSEC_FOUND) OR (PARSEC_GIVEN_BY_USER) )

  # Looking for include
  # -------------------

  # Add system include paths to search include
  # ------------------------------------------
  unset(_inc_env)
  set(ENV_PARSEC_DIR "$ENV{PARSEC_DIR}")
  set(ENV_PARSEC_INCDIR "$ENV{PARSEC_INCDIR}")
  if(ENV_PARSEC_INCDIR)
    list(APPEND _inc_env "${ENV_PARSEC_INCDIR}")
  elseif(ENV_PARSEC_DIR)
    list(APPEND _inc_env "${ENV_PARSEC_DIR}")
    list(APPEND _inc_env "${ENV_PARSEC_DIR}/include")
    list(APPEND _inc_env "${ENV_PARSEC_DIR}/include/parsec")
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

  # Try to find the parsec headers in the given paths
  # -------------------------------------------------

  # create list of headers to find
  set(PARSEC_hdrs_to_find "parsec_config.h" "parsec.h")

  # call cmake macro to find the header path
  if(PARSEC_INCDIR)
    foreach(parsec_hdr ${PARSEC_hdrs_to_find})
      set(PARSEC_${parsec_hdr}_INCLUDE_DIRS "PARSEC_${parsec_hdr}_INCLUDE_DIRS-NOTFOUND")
      find_path(PARSEC_${parsec_hdr}_INCLUDE_DIRS
        NAMES ${parsec_hdr}
        HINTS ${PARSEC_INCDIR})
    endforeach()
  else()
    if(PARSEC_DIR)
      set(PARSEC_${parsec_hdr}_INCLUDE_DIRS "PARSEC_${parsec_hdr}_INCLUDE_DIRS-NOTFOUND")
      foreach(parsec_hdr ${PARSEC_hdrs_to_find})
        find_path(PARSEC_${parsec_hdr}_INCLUDE_DIRS
          NAMES ${parsec_hdr}
          HINTS ${PARSEC_DIR}
          PATH_SUFFIXES "include")
      endforeach()
    else()
      foreach(parsec_hdr ${PARSEC_hdrs_to_find})
        set(PARSEC_${parsec_hdr}_INCLUDE_DIRS "PARSEC_${parsec_hdr}_INCLUDE_DIRS-NOTFOUND")
        find_path(PARSEC_${parsec_hdr}_INCLUDE_DIRS
          NAMES ${parsec_hdr}
          HINTS ${_inc_env})
      endforeach()
    endif()
  endif()

  # If found, add path to cmake variable
  # ------------------------------------
  set(PARSEC_INCLUDE_DIRS "")
  foreach(parsec_hdr ${PARSEC_hdrs_to_find})

    if (PARSEC_${parsec_hdr}_INCLUDE_DIRS)
      # set cmake variables using the pkg-config naming convention
      list(APPEND PARSEC_INCLUDE_DIRS "${PARSEC_${parsec_hdr}_INCLUDE_DIRS}" )
    else ()
      if(NOT PARSEC_FIND_QUIETLY)
        message(STATUS "Looking for parsec -- ${parsec_hdr} not found")
      endif()
    endif ()
    mark_as_advanced(PARSEC_${parsec_hdr}_INCLUDE_DIRS)

  endforeach(parsec_hdr ${PARSEC_hdrs_to_find})

  if (PARSEC_INCLUDE_DIRS)
    list(REMOVE_DUPLICATES PARSEC_INCLUDE_DIRS)
  endif ()

  # Looking for lib
  # ---------------

  set(PARSEC_LIBRARIES "")
  set(PARSEC_LIBRARY_DIRS "")

  # Add system library paths to search lib
  # --------------------------------------
  unset(_lib_env)
  set(ENV_PARSEC_LIBDIR "$ENV{PARSEC_LIBDIR}")
  if(ENV_PARSEC_LIBDIR)
    list(APPEND _lib_env "${ENV_PARSEC_LIBDIR}")
  elseif(ENV_PARSEC_DIR)
    list(APPEND _lib_env "${ENV_PARSEC_DIR}")
    list(APPEND _lib_env "${ENV_PARSEC_DIR}/lib")
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

  # Try to find the parsec libs in the given paths
  # ----------------------------------------------

  # create list of libs to find
  set(PARSEC_libs_to_find "parsec" "parsec-base" "parsec_distribution" "parsec_distribution_matrix")

  # call cmake macro to find the lib path
  if(PARSEC_LIBDIR)
    foreach(parsec_lib ${PARSEC_libs_to_find})
      set(PARSEC_${parsec_lib}_LIBRARY "PARSEC_${parsec_lib}_LIBRARY-NOTFOUND")
      find_library(PARSEC_${parsec_lib}_LIBRARY
        NAMES ${parsec_lib}
        HINTS ${PARSEC_LIBDIR})
    endforeach()
  else()
    if(PARSEC_DIR)
      foreach(parsec_lib ${PARSEC_libs_to_find})
        set(PARSEC_${parsec_lib}_LIBRARY "PARSEC_${parsec_lib}_LIBRARY-NOTFOUND")
        find_library(PARSEC_${parsec_lib}_LIBRARY
          NAMES ${parsec_lib}
          HINTS ${PARSEC_DIR}
          PATH_SUFFIXES lib lib32 lib64)
      endforeach()
    else()
      foreach(parsec_lib ${PARSEC_libs_to_find})
        set(PARSEC_${parsec_lib}_LIBRARY "PARSEC_${parsec_lib}_LIBRARY-NOTFOUND")
        find_library(PARSEC_${parsec_lib}_LIBRARY
          NAMES ${parsec_lib}
          HINTS ${_lib_env})
      endforeach()
    endif()
  endif()

  # If found, add path to cmake variable
  # ------------------------------------
  foreach(parsec_lib ${PARSEC_libs_to_find})

    if (PARSEC_${parsec_lib}_LIBRARY)

      get_filename_component(${parsec_lib}_lib_path ${PARSEC_${parsec_lib}_LIBRARY} PATH)
      # set cmake variables (respects naming convention)
      list(APPEND PARSEC_LIBRARIES "${PARSEC_${parsec_lib}_LIBRARY}")
      list(APPEND PARSEC_LIBRARY_DIRS "${${parsec_lib}_lib_path}")

    else (PARSEC_${parsec_lib}_LIBRARY)

      if(NOT PARSEC_FIND_QUIETLY)
        message(STATUS "Looking for parsec -- lib ${parsec_lib} not found")
      endif()

    endif (PARSEC_${parsec_lib}_LIBRARY)

    mark_as_advanced(PARSEC_${parsec_lib}_LIBRARY)

  endforeach(parsec_lib ${PARSEC_libs_to_find})

  list(REMOVE_DUPLICATES PARSEC_LIBRARIES)
  if (PARSEC_LIBRARY_DIRS)
    list(REMOVE_DUPLICATES PARSEC_LIBRARY_DIRS)
  endif ()

  # Looking for parsec compilers
  # ----------------------------

  # Add system bin paths to search drivers
  # --------------------------------------
  unset(_bin_env)
  if(ENV_PARSEC_DIR)
    list(APPEND _bin_env "${ENV_PARSEC_DIR}")
    list(APPEND _bin_env "${ENV_PARSEC_DIR}/bin")
  else()
    string(REPLACE ":" ";" _bin_env "$ENV{PATH}")
  endif()
  list(REMOVE_DUPLICATES _bin_env)

  # create list of binaries to find
  set(PARSEC_bins_to_find "parsec_ptgpp")

  # call cmake macro to find the bin path
  if(PARSEC_DIR)
    foreach(parsec_bin ${PARSEC_bins_to_find})
      set(PARSEC_${parsec_bin}_BIN_DIR "PARSEC_${parsec_bin}_BIN_DIR-NOTFOUND")
      find_path(PARSEC_${parsec_bin}_BIN_DIR
        NAMES ${parsec_bin}
        HINTS ${PARSEC_DIR}
        PATH_SUFFIXES "bin")
    endforeach()
  else()
    foreach(parsec_bin ${PARSEC_bins_to_find})
      set(PARSEC_${parsec_bin}_BIN_DIR "PARSEC_${parsec_bin}_BIN_DIR-NOTFOUND")
      find_path(PARSEC_${parsec_bin}_BIN_DIR
        NAMES ${parsec_bin}
        HINTS ${_bin_env})
    endforeach()
  endif()
  if (PARSEC_parsec_ptgpp_BIN_DIR)
    if (NOT PARSEC_FIND_QUIETLY)
      message(STATUS "Look for PARSEC - compiler parsec_ptgpp found in ${PARSEC_parsec_ptgpp_BIN_DIR}")
    endif()
  else()
    if (PARSEC_FIND_REQUIRED)
      message(FATAL_ERROR "Look for PARSEC - compiler parsec_ptgpp not found while required")
    endif()
  endif()

endif( (NOT PKG_CONFIG_EXECUTABLE) OR (PKG_CONFIG_EXECUTABLE AND NOT PARSEC_FOUND) OR (PARSEC_GIVEN_BY_USER) )

# check a function to validate the find
if(PARSEC_LIBRARIES)

  set(REQUIRED_FLAGS)
  set(REQUIRED_LDFLAGS)
  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # PARSEC
  if (PARSEC_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${PARSEC_INCLUDE_DIRS}")
  endif()
  if (PARSEC_CFLAGS_OTHER)
    set(REQUIRED_FLAGS "${PARSEC_CFLAGS_OTHER}")
  endif()
  if (PARSEC_LDFLAGS_OTHER)
    set(REQUIRED_LDFLAGS "${PARSEC_LDFLAGS_OTHER}")
  endif()
  foreach(libdir ${PARSEC_LIBRARY_DIRS})
    if (libdir)
      list(APPEND REQUIRED_LIBDIRS "${libdir}")
    endif()
  endforeach()
  set(REQUIRED_LIBS "${PARSEC_LIBRARIES}")
  # HWLOC
  if (HWLOC_FOUND AND PARSEC_LOOK_FOR_HWLOC)
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
  # MPI
  if (MPI_FOUND AND PARSEC_LOOK_FOR_MPI)
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
  # CUDA
  if (CUDA_FOUND AND PARSEC_LOOK_FOR_CUDA)
    if (CUDA_INCLUDE_DIRS)
      list(APPEND REQUIRED_INCDIRS "${CUDA_INCLUDE_DIRS}")
    endif()
    if (CUDA_LIBRARY_DIRS)
      list(APPEND REQUIRED_LIBDIRS "${CUDA_LIBRARY_DIRS}")
    endif()
    list(APPEND REQUIRED_LIBS "${CUDA_CUBLAS_LIBRARIES};${CUDA_CUDART_LIBRARY};${CUDA_CUDA_LIBRARY}")
  endif()
  # EXTRA LIBS such that pthread, m, rt, dl
  list(APPEND REQUIRED_LIBS ${PARSEC_EXTRA_LIBRARIES})

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
  unset(PARSEC_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(parsec_init PARSEC_WORKS)
  mark_as_advanced(PARSEC_WORKS)

  if(PARSEC_WORKS)
    # save link with dependencies
    set(PARSEC_LIBRARIES_DEP "${REQUIRED_LIBS}")
    set(PARSEC_LIBRARY_DIRS_DEP "${REQUIRED_LIBDIRS}")
    set(PARSEC_INCLUDE_DIRS_DEP "${REQUIRED_INCDIRS}")
    set(PARSEC_CFLAGS_OTHER_DEP "${REQUIRED_FLAGS}")
    set(PARSEC_LDFLAGS_OTHER_DEP "${REQUIRED_LDFLAGS}")
  else()
    if(NOT PARSEC_FIND_QUIETLY)
      message(STATUS "Looking for parsec : test of parsec_init fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "CMAKE_REQUIRED_FLAGS: ${CMAKE_REQUIRED_FLAGS}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
      message(STATUS "Maybe PARSEC is linked with specific libraries. "
        "Have you tried with COMPONENTS (HWLOC, CUDA, MPI)? "
        "See the explanation in FindPARSEC.cmake.")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(PARSEC_LIBRARIES)

if (PARSEC_LIBRARIES)
  if (PARSEC_LIBRARY_DIRS)
    foreach(dir ${PARSEC_LIBRARY_DIRS})
      if ("${dir}" MATCHES "parsec")
        set(first_lib_path "${dir}")
      endif()
    endforeach()
  else()
    list(GET PARSEC_LIBRARIES 0 first_lib)
    get_filename_component(first_lib_path "${first_lib}" PATH)
    set(PARSEC_LIBRARY_DIRS "${first_lib_path}")
  endif()
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(PARSEC_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of PARSEC library" FORCE)
  else()
    set(PARSEC_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of PARSEC library" FORCE)
  endif()
endif()
mark_as_advanced(PARSEC_DIR)
mark_as_advanced(PARSEC_DIR_FOUND)

# check that PARSEC has been found
# --------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARSEC DEFAULT_MSG
  PARSEC_LIBRARIES
  PARSEC_parsec_ptgpp_BIN_DIR
  PARSEC_WORKS)

if ( PARSEC_parsec_ptgpp_BIN_DIR )
  find_program(PARSEC_PARSEC_PTGPP
    NAMES parsec_ptgpp
    HINTS ${PARSEC_parsec_ptgpp_BIN_DIR})
else()
  set(PARSEC_PARSEC_PTGPP "PARSEC_PARSEC_PTGPP-NOTFOUND")
endif()
