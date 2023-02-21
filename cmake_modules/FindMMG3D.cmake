## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Inria - IMB (Université de Bordeaux) - LJLL (UPMC), 2004- .
##
##  mmg is free software: you can redistribute it and/or modify it
##  under the terms of the GNU Lesser General Public License as published
##  by the Free Software Foundation, either version 3 of the License, or
##  (at your option) any later version.
##
##  mmg is distributed in the hope that it will be useful, but WITHOUT
##  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
##  License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License and of the GNU General Public License along with mmg (in
##  files COPYING.LESSER and COPYING). If not, see
##  <http://www.gnu.org/licenses/>. Please read their terms carefully and
##  use this copy of the mmg distribution only if you accept them.
## =============================================================================
#
# - Find mmg3d include dirs and libraries
# Use this module by invoking find_package with the form:
#  find_package(MMG3D
#              [REQUIRED]             # Fail with error if mmg3d not found
#              )
#
# This module finds headers and mmg libraries.
# Results are reported in variables:
#  MMG3D_FOUND           - True if headers and requested libraries were found
#  MMG3D_INCLUDE_DIRS    - mmg3d include directories
#  MMG3D_LIBRARY_DIRS    - Link directories for mmg libraries
#  MMG3D_LIBRARIES       - mmg component libraries to be linked
#
# The user can give specific paths where to find the libraries adding cmake
# options at configure (ex: cmake path/to/project -DMMG_DIR=path/to/mmg/install):
#  MMG_DIR                  - Where to find the base directory of mmg
#  MMG_INCDIR               - Path toward Mmg installation (should at least contains mmg/common filder)
#  MMG_LIBDIR               - Where to find the library files
# The module can also look for the following environment variables if paths
# are not given as cmake variable: MMG_DIR, MMG_INCDIR, MMG_LIBDIR

if (NOT MMG3D_FOUND)
  set(MMG_DIR "" CACHE PATH "Installation directory of MMG library")
  if (NOT MMG_FIND_QUIETLY)
    message(STATUS "A cache variable, namely MMG_DIR, has been set to specify the install directory of MMG")
  endif()
endif()

# Looking for include
# -------------------

# Add system include paths to search include
# ------------------------------------------
unset(_inc_env)
set(ENV_MMG_DIR "$ENV{MMG_DIR}")
set(ENV_MMG_INCDIR "$ENV{MMG_INCDIR}")

if(ENV_MMG_INCDIR)
  list(APPEND _inc_env "${ENV_MMG_INCDIR}")
elseif(ENV_MMG_DIR)
  list(APPEND _inc_env "${ENV_MMG_DIR}")
  list(APPEND _inc_env "${ENV_MMG_DIR}/include")
  list(APPEND _inc_env "${ENV_MMG_DIR}/include/mmg")
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
list(APPEND _inc_env "${CMAKE_PLATFORM_IMPLICIT_INCLUDE_DIRECTORIES}")
list(APPEND _inc_env "${CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES}")
list(REMOVE_DUPLICATES _inc_env)


# Try to find the header in the given paths
# -------------------------------------------------
# call cmake macro to find the header path
if(MMG_INCDIR)
  set(MMG3D_libmmgtypes.h_DIRS "MMG3D_libmmgtypes.h_DIRS-NOTFOUND")
  find_path(MMG3D_libmmgtypes.h_DIRS
    NAMES libmmgtypes.h
    HINTS ${MMG_INCDIR}
    PATH_SUFFIXES "mmg/mmg3d" "mmg3d" "mmg/common")
else()
  if(MMG_DIR)
    set(MMG3D_libmmgtypes.h_DIRS "MMG3D_libmmgtypes.h_DIRS-NOTFOUND")
    find_path(MMG3D_libmmgtypes.h_DIRS
      NAMES libmmgtypes.h
      HINTS ${MMG_DIR}/include
      PATH_SUFFIXES "mmg" "mmg/common")
  else()
    set(MMG3D_libmmgtypes.h_DIRS "MMG3D_libmmgtypes.h_DIRS-NOTFOUND")
    find_path(MMG3D_libmmgtypes.h_DIRS
      NAMES libmmgtypes.h
      HINTS ${_inc_env}
      PATH_SUFFIXES  "mmg" "mmg/common")
  endif()
endif()
STRING(REGEX REPLACE "(mmg/mmg3d)|(mmg/common)" "" MMG3D_libmmgtypes.h_DIRS "${MMG3D_libmmgtypes.h_DIRS}" )

mark_as_advanced(MMG3D_libmmgtypes.h_DIRS)


# If found, add path to cmake variable
# ------------------------------------
if (MMG3D_libmmgtypes.h_DIRS AND NOT MMG_FIND_QUIETLY)
  set(MMG3D_INCLUDE_DIRS "${MMG3D_libmmgtypes.h_DIRS}")
else ()
  set(MMG3D_INCLUDE_DIRS "MMG3D_INCLUDE_DIRS-NOTFOUND")
  if(NOT MMG_FIND_QUIETLY)
    message(STATUS "Looking for mmg3d -- libmmgtypes.h not found")
  endif()
endif()


# Looking for lib
# ---------------

# Add system library paths to search lib
# --------------------------------------
unset(_lib_env)
set(ENV_MMG_LIBDIR "$ENV{MMG_LIBDIR}")
if(ENV_MMG_LIBDIR)
  list(APPEND _lib_env "${ENV_MMG_LIBDIR}")
elseif(ENV_MMG_DIR)
  list(APPEND _lib_env "${ENV_MMG_DIR}")
  list(APPEND _lib_env "${ENV_MMG_DIR}/lib")
else()
  if(WIN32)
    string(REPLACE ":" ";" _lib_env "$ENV{LIB}")
  else()
    if(APPLE)
      string(REPLACE ":" ";" _lib_env "$ENV{DYLD_LIBRARY_PATH}")
    else()
      string(REPLACE ":" ";" _lib_env "$ENV{LD_LIBRARY_PATH}")
    endif()
    list(APPEND _lib_env "${CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES}")
    list(APPEND _lib_env "${CMAKE_C_IMPLICIT_LINK_DIRECTORIES}")
  endif()
endif()
list(REMOVE_DUPLICATES _lib_env)

# Try to find the mmg3d lib in the given paths
# ----------------------------------------------
# call cmake macro to find the lib path
if(MMG_LIBDIR)
  set(MMG3D_mmg3d_LIBRARY "MMG3D_mmg3d_LIBRARY-NOTFOUND")
  find_library(MMG3D_mmg3d_LIBRARY
    NAMES mmg3d
    HINTS ${MMG_LIBDIR})
else()
  if(MMG_DIR)
    set(MMG3D_mmg3d_LIBRARY "MMG3D_mmg3d_LIBRARY-NOTFOUND")

    find_library(MMG3D_mmg3d_LIBRARY
      NAMES mmg3d
      HINTS ${MMG_DIR}
      PATH_SUFFIXES lib lib32 lib64)
  else()
    set(MMG3D_mmg3d_LIBRARY "MMG3D_mmg3d_LIBRARY-NOTFOUND")
    find_library(MMG3D_mmg3d_LIBRARY
      NAMES mmg3d
      HINTS ${_lib_env})
  endif()
endif()
mark_as_advanced(MMG3D_mmg3d_LIBRARY)

# If found, add path to cmake variable
# ------------------------------------
if (MMG3D_mmg3d_LIBRARY)
  get_filename_component(mmg_lib_path "${MMG3D_mmg3d_LIBRARY}" PATH)
  # set cmake variables
  set(MMG3D_LIBRARIES    "${MMG3D_mmg3d_LIBRARY}")
  set(MMG3D_LIBRARY_DIRS "${mmg_lib_path}")
else ()
  set(MMG3D_LIBRARIES    "MMG3D_LIBRARIES-NOTFOUND")
  set(MMG3D_LIBRARY_DIRS "MMG3D_LIBRARY_DIRS-NOTFOUND")
  if (NOT MMG_FIND_QUIETLY)
    message(STATUS "Looking for mmg3d -- lib mmg3d not found")
  endif()
endif ()

# check a function to validate the find
if(MMG3D_LIBRARIES)

  set(REQUIRED_INCDIRS)
  set(REQUIRED_LIBDIRS)
  set(REQUIRED_LIBS)

  # MMG3D
  if (MMG3D_INCLUDE_DIRS)
    set(REQUIRED_INCDIRS "${MMG3D_INCLUDE_DIRS}")
  endif()
  if (MMG3D_LIBRARY_DIRS)
    set(REQUIRED_LIBDIRS "${MMG3D_LIBRARY_DIRS}")
  endif()
  set(REQUIRED_LIBS "${MMG3D_LIBRARIES}")
  # m
  find_library(M_LIBRARY NAMES m)
  mark_as_advanced(M_LIBRARY)
  if(M_LIBRARY)
    list(APPEND REQUIRED_LIBS "-lm")
  endif()

  # set required libraries for link
  set(CMAKE_REQUIRED_INCLUDES "${REQUIRED_INCDIRS}")
  set(CMAKE_REQUIRED_LIBRARIES)
  foreach(lib_dir ${REQUIRED_LIBDIRS})
    list(APPEND CMAKE_REQUIRED_LIBRARIES "-L${lib_dir}")
  endforeach()
  list(APPEND CMAKE_REQUIRED_LIBRARIES "${REQUIRED_LIBS}")
  string(REGEX REPLACE "^ -" "-" CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES}")

  # test link
  unset(MMG3D_WORKS CACHE)
  include(CheckFunctionExists)
  check_function_exists(MMG5_loadMshMesh_part1 MMG3D_WORKS)
  mark_as_advanced(MMG3D_WORKS)

  if(NOT MMG3D_WORKS)
    if(NOT MMG_FIND_QUIETLY)
      message(STATUS "Looking for MMG3D : test of MMG5_loadMshMesh_part1 with MMG3D library fails")
      message(STATUS "CMAKE_REQUIRED_LIBRARIES: ${CMAKE_REQUIRED_LIBRARIES}")
      message(STATUS "CMAKE_REQUIRED_INCLUDES: ${CMAKE_REQUIRED_INCLUDES}")
      message(STATUS "Check in CMakeFiles/CMakeError.log to figure out why it fails")
    endif()
  endif()
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LIBRARIES)
endif(MMG3D_LIBRARIES)

if (MMG3D_LIBRARIES)
  list(GET MMG3D_LIBRARIES 0 first_lib)
  get_filename_component(first_lib_path "${first_lib}" PATH)
  if (${first_lib_path} MATCHES "/lib(32|64)?$")
    string(REGEX REPLACE "/lib(32|64)?$" "" not_cached_dir "${first_lib_path}")
    set(MMG_DIR_FOUND "${not_cached_dir}" CACHE PATH "Installation directory of MMG library" FORCE)
  else()
    set(MMG_DIR_FOUND "${first_lib_path}" CACHE PATH "Installation directory of MMG library" FORCE)
  endif()
endif()
mark_as_advanced(MMG_DIR_FOUND)

# check that MMG has been found
# ----------------------------------
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MMG3D DEFAULT_MSG
  MMG3D_LIBRARIES
  MMG3D_WORKS)