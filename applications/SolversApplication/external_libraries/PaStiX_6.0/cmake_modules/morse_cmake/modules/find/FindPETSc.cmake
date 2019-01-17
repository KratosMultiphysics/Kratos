#
# @copyright (c) 2018 Inria. All rights reserved.
#
# Marc Fuentes
# Florent Pruvost
#
# https://gitlab.inria.fr/sed-bso/findPetsc
#
# INPUT
# -----
#
# PETSc is not installed in a standard way on Unix systems so that
# this module requires hints to know where PETSc is installed. Please
# give the installation directory (contains ./include/petsc.h, ./lib/, etc):
# 1. by setting the PETSC_DIR variable
#    a. as an environment variable, e.g.
#       $ export PETSC_DIR=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real
#    b. or as an CMake variable, e.g.
#       $ cmake .. -DPETSC_DIR==/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real
# 2. or by using the pkg-config mechanism, e.g.
#    $ export PKG_CONFIG_PATH=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real/lib/pkgconfig:$PKG_CONFIG_PATH
#
# OUTPUT
# -------
# PETSC_INCLUDE_DIRS - the PETSc include directories
# PETSC_LIBRARIES    - Link these to use PETSc
# PETSC_LIBRARY_DIRS - Link these to use PETSc
# PETSC_MPIEXEC - Executable for running MPI programs
#
# if pkg-config is used i.e. pkgconfig installed, PETSc.pc file path
# in the PKG_CONFIG_PATH environment variable and PETSC_DIR not set
# then the following variables are set (or empty)
#
# PETSC_VERSION           ... version of the module
# PETSC_PREFIX            ... prefix-directory of the module
# PETSC_INCLUDEDIR        ... include-dir of the module
# PETSC_LIBDIR            ... lib-dir of the module
#
# <PREFIX>_FOUND          ... set to 1 if module(s) exist
# <PREFIX>_LIBRARIES      ... only the libraries (w/o the '-l')
# <PREFIX>_LIBRARY_DIRS   ... the paths of the libraries (w/o the '-L')
# <PREFIX>_LDFLAGS        ... all required linker flags
# <PREFIX>_LDFLAGS_OTHER  ... all other linker flags
# <PREFIX>_INCLUDE_DIRS   ... the '-I' preprocessor flags (w/o the '-I')
# <PREFIX>_CFLAGS         ... all required cflags
# <PREFIX>_CFLAGS_OTHER   ... the other compiler flags
#
# <PREFIX> = PETSC        for common case
# <PREFIX> = PETSC_STATIC for static linking
#
# PETSC_FOUND_WITH_PKGCONFIG - True if found with pkg-config
#
# find_package(PETSc [QUIET] [REQUIRED])
#
# Setting these changes the behavior of the search
# PETSC_DIR - directory in which PETSc is installed
# PETSC_ARCH - build architecture

# create a cmake cache variable
set(PETSC_DIR "" CACHE PATH "Installation directory of PETSC library")
if (NOT PETSc_FIND_QUIETLY AND NOT PETSC_DIR)
  message(STATUS "A cache variable, namely PETSC_DIR, has been set
  to specify a custom installation directory of PETSC")
endif()

# Use pkg-config to detect include/library dirs (if pkg-config is available)
# -------------------------------------------------------------------------------------
include(FindPkgConfig)
find_package(PkgConfig QUIET)
if( PKG_CONFIG_EXECUTABLE AND NOT PETSC_DIR )
  pkg_search_module(PETSC PETSc)
  if (NOT PETSc_FIND_QUIETLY)
    if (PETSC_FOUND AND PETSC_LIBRARIES)
      message(STATUS "Looking for PETSC - found using PkgConfig")
    else()
      message(STATUS "Looking for PETSC - not found using PkgConfig."
        "\n   Perhaps you should add the directory containing PETSc.pc to"
        "\n   the PKG_CONFIG_PATH environment variable.")
    endif()
  endif()
  set(PETSC_DIR "${PETSC_PREFIX}")
  if (PETSC_FOUND AND PETSC_LIBRARIES)
    set(PETSC_FOUND_WITH_PKGCONFIG "TRUE")
    find_pkgconfig_libraries_absolute_path(PETSC)
  else()
    set(PETSC_FOUND_WITH_PKGCONFIG "FALSE")
  endif()
endif()

# consider using the env. var. PETSC_DIR if not directly given through the CMake cache var.
if (NOT PETSC_DIR AND DEFINED ENV{PETSC_DIR})
  set(PETSC_DIR "$ENV{PETSC_DIR}")
endif()

if (PETSC_DIR)
    if (EXISTS ${PETSC_DIR})
      if (EXISTS ${PETSC_DIR}/include/petsc.h)
          if (NOT PETSc_FIND_QUIETLY)
            message(STATUS "PETSC_DIR = ${PETSC_DIR} contains include/petsc.h")
          endif()
        else()
          if (PETSc_FIND_REQUIRED)
            message(FATAL_ERROR "include/petsc.h not found in PETSC_DIR = ${PETSC_DIR}")
          endif()
        endif()
      else()
        if (PETSc_FIND_REQUIRED)
          message(FATAL_ERROR "PETSC_DIR defined, but ${PETSC_DIR} does not exist")
        endif()
    endif()
else()
  if (PETSc_FIND_REQUIRED)
    message(FATAL_ERROR "\
PETSc is not installed in a standard way on Unix systems so that
this module requires hints to know where PETSc is installed. Please
give the installation directory (contains ./include/petsc.h, ./lib/, etc):
1. by setting the PETSC_DIR variable
   a. as an environment variable, e.g.
      $ export PETSC_DIR=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real
   b. or as an CMake variable, e.g.
      $ cmake .. -DPETSC_DIR==/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real
2. or by using the pkg-config mechanism, e.g.
   $ export PKG_CONFIG_PATH=/usr/lib/petscdir/3.6.2/x86_64-linux-gnu-real/lib/pkgconfig:$PKG_CONFIG_PATH\
    ")
  endif()
endif()

find_file(petscconf NAMES petscconf.h HINTS ${PETSC_DIR}/include)
if (petscconf)
  if (NOT PETSc_FIND_QUIETLY)
    message(STATUS "petscconf.h is located ad ${petscconf}")
  endif()
else()
  if (PETSc_FIND_REQUIRED)
    message(FATAL_ERROR "\petscconf.h not found
    we must find it into PETSC_DIR/include/")
  endif()
endif()

if (NOT PETSC_ARCH)
  if (petscconf)
     file(READ ${petscconf} contents)
     string(REGEX MATCHALL "#define *PETSC_ARCH *\".+\" *\n" foundArch ${contents})
     string(REGEX REPLACE "#define *PETSC_ARCH *\"(.+)\" *\n" "\\1" archFromFile ${foundArch})
     if (NOT PETSc_FIND_QUIETLY)
       message(STATUS "PETSC_ARCH taken from petscconf.h = ${archFromFile}")
     endif()
     if(archFromFile)
         set(PETSC_ARCH ${archFromFile})
     endif(archFromFile)
  endif(petscconf)
endif (NOT PETSC_ARCH)
if (NOT PETSc_FIND_QUIETLY)
  message(STATUS "PETSC_ARCH  = ${PETSC_ARCH}")
endif()

find_file(petscvariables NAMES petscvariables HINTS ${PETSC_DIR}/lib/petsc/conf/ ${PETSC_DIR}/conf)
if (petscvariables)
  if (NOT PETSc_FIND_QUIETLY)
    message(STATUS "petscvariables = ${petscvariables}")
  endif()
else()
  if (PETSc_FIND_REQUIRED)
    message(FATAL_ERROR "petscvariables not found")
  endif()
endif()

# function extracting the value of variable defined in a file under the form
# NAME_VAR = VARIABLE_CONTENT
# note that VARIABLE_CONTENT is assumed not to contain ";"
function(get_variable var filename out)
    if (NOT EXISTS ${filename})
        message(FATAL_ERROR "${filename} does not exist")
    endif()
    file(READ ${filename} contents)
    string(REGEX REPLACE "\n" ";" contents "${contents}") # we divide line by lines because it did work without doing that
    foreach(Line ${contents})
        string(REGEX MATCHALL "^${var} = .*" foundVar ${Line})
        if (foundVar)
            string(REGEX REPLACE "^${var} = (.+)" "\\1" value ${foundVar})
            break()
        endif()
    endforeach()
    if(foundVar)
      set(${out} "${value}" PARENT_SCOPE)
    else()
        message(FATAL_ERROR "${var} not found in ${filename}")
    endif()
endfunction(get_variable)

get_variable(PETSC_WITH_EXTERNAL_LIB ${petscvariables} PETSC_LIBRARIES_RAW)
get_variable(MPIEXEC ${petscvariables} PETSC_MPIEXEC)
get_variable(PETSC_CC_INCLUDES ${petscvariables}  PETSC_INCLUDES_RAW)

# extract include paths from petscvariables
string(REGEX REPLACE "-I" "" PETSC_INCLUDES_RAW ${PETSC_INCLUDES_RAW} )
string(REGEX REPLACE " " ";" PETSC_INCLUDE_DIRS ${PETSC_INCLUDES_RAW} )

# extract libraries  from petscvariables
string(REGEX MATCHALL " -l([^ ]+) " foundLibs ${PETSC_LIBRARIES_RAW})
string(REGEX REPLACE " -l" "" PETSC_LIBRARIES ${foundLibs})
string(REGEX REPLACE " " ";" PETSC_LIBRARIES ${PETSC_LIBRARIES})

# extract library directories  from petscvariables
string(REGEX MATCHALL " -L([^ ]+) " foundLibPaths ${PETSC_LIBRARIES_RAW})
string(REGEX REPLACE " -L" "" PETSC_LIBRARY_DIRS ${foundLibPaths})
string(REGEX REPLACE " " ";" PETSC_LIBRARY_DIRS ${PETSC_LIBRARY_DIRS})

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PETSc
  "PETSc could not be found. Be sure to set PETSC_DIR."
  PETSC_MPIEXEC PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_LIBRARY_DIRS)
