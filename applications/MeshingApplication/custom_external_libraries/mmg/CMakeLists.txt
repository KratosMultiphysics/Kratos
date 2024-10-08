## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.12)

# if PROJECT_NAME is defined, mmg is a subproject
if(DEFINED PROJECT_NAME)
  set(MMG_SUBPROJECT ON)
endif()

INCLUDE(CMakeDependentOption)

PROJECT (mmg)

LIST(APPEND CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake/modules
  ${PROJECT_SOURCE_DIR}/cmake/testing )

INCLUDE(macros)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(MMGCOMMON_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src/common)
SET(MMGCOMMON_BINARY_DIR ${PROJECT_BINARY_DIR}/src/common)

FILE(MAKE_DIRECTORY ${MMGCOMMON_BINARY_DIR})

# Path toward the continuous integration directory
SET ( CI_DIR ${PROJECT_SOURCE_DIR}/ci_tests
  CACHE PATH "Path toward the continuous integration directory" )
MARK_AS_ADVANCED ( CI_DIR )


include(CMakePackageConfigHelpers)
# various variables to help standard install directories
include(GNUInstallDirs)


###############################################################################
#####
#####         Release version and date
#####
###############################################################################

SET (CMAKE_RELEASE_VERSION_MAJOR "5" )
SET (CMAKE_RELEASE_VERSION_MINOR "7" )
SET (CMAKE_RELEASE_VERSION_PATCH "0" )
SET (CMAKE_RELEASE_DATE "Dec. 13, 2022" )

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

###############################################################################
#####
#####  Set preprocessor flags and manage automatically generated headers data
#####
###############################################################################
# Architecture detection
IF (NOT WIN32)
# Set preprocessor flags to say that we are posix and gnu compatible
  SET (DEF_POSIX "#define MMG_POSIX")
  SET (DEF_GNU "#define MMG_GNU")
ELSEIF(MINGW)
  SET (DEF_GNU "#define MMG_GNU")
ENDIF ()

############################################################################
#####
#####         Integer type
#####
############################################################################

# type of integer
SET (MMG5_INT "int32_t" CACHE STRING "Integer type")
IF (${MMG5_INT} STREQUAL "int64_t")
  SET ( DEF_MMG5_INT "#define MMG5_int int64_t")
  SET ( DEF_MMG5_PRId "#define MMG5_PRId PRId64")
  SET ( DEF_MMG5F_INT "#define MMG5F_INT SELECTED_INT_KIND(18)")
  SET ( DEF_MMG_SWPBIN "#define MMG5_SWAPBIN(x) MMG5_swapbin_int(x)")
  SET ( DEF_MMG_ABS "#define MMG5_abs(x) (int64_t)labs((int64_t)x)")
ELSE ()
  SET ( DEF_MMG5_INT "#define MMG5_int int32_t")
  SET ( DEF_MMG5_PRId "#define MMG5_PRId PRId32")
  SET ( DEF_MMG5F_INT "#define MMG5F_INT SELECTED_INT_KIND(9)")
  SET ( DEF_MMG_SWPBIN "#define MMG5_SWAPBIN(x) MMG5_swapbin(x)")
  SET ( DEF_MMG_ABS "#define MMG5_abs(x) abs(x)")
ENDIF ()

# Create mmgcmakedefines.h file with the architecture infos and the definition
# ot the USE_POINTMAP flag. Create also mmgcmakedefinesf.h file for automatic
# definition of integer size for Fortran API
###############################################################################
SET(MMG_DYN_LIB ${BUILD_SHARED_LIBS} CACHE INTERNAL "" FORCE)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgcmakedefines.h.in
  ${PROJECT_BINARY_DIR}/src/common/mmgcmakedefines.h @ONLY)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgcmakedefinesf.h.in
  ${PROJECT_BINARY_DIR}/src/common/mmgcmakedefinesf.h @ONLY)

# Create mmgversion.h file with the release infos
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgversion.h.in
  ${PROJECT_BINARY_DIR}/src/common/mmgversion.h @ONLY)

# Print git hash to source file
IF (NOT WIN32 OR MINGW)

  ADD_CUSTOM_TARGET(GenerateGitHashMmg
    COMMAND ./git_log_mmg.sh ${PROJECT_SOURCE_DIR} ${MMGCOMMON_BINARY_DIR}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/scripts/
    COMMENT "Getting git commit hash"
    BYPRODUCTS ${MMGCOMMON_BINARY_DIR}/git_log_mmg.h
  )
ENDIF ()

###############################################################################
#####
#####         CMake flags / variables
#####
###############################################################################
# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT)

include(add_build_types)

# List CMake build types
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  SET(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo" "Maintainer" "RelWithAssert")
ELSEIF( CMAKE_CONFIGURATION_TYPES )
  LIST(APPEND CMAKE_CONFIGURATION_TYPES Maintainer RelWithAssert )
ENDIF()

# Explicitly set the DNDEBUG flag in case the user or a parent project overrides
# it.
if (NOT CMAKE_BUILD_TYPE MATCHES Debug AND NOT CMAKE_BUILD_TYPE MATCHES RelWithAssert)
    add_definitions(-DNDEBUG)
endif()

# Do not print warnings about the use of secure function with MSVC and the use
# of POSIX strdup
if ( MSVC )
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-Dstrdup=_strdup)
endif()

###############################################################################
#####
#####         Choose executable target to compile
#####
###############################################################################
SET ( BUILD "MMG" CACHE STRING
  "\"MMG3D\", \"MMGS\",\"MMG2D\", or \"MMG\" to compile both MMGS, MMG3D and MMG2D" )

IF ( BUILD MATCHES "MMGS" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       ON  )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG3D" )
  SET(BUILD_MMG3D      ON  )
  SET(BUILD_MMGS       OFF )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG2D" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       OFF  )
  SET(BUILD_MMG2D      ON  )
ELSEIF ( BUILD MATCHES "MMG" )
  SET(BUILD_MMG        ON  )
  SET(BUILD_MMGS       ON )
  SET(BUILD_MMG3D      ON )
  SET(BUILD_MMG2D      ON  )
ELSE ( )
  MESSAGE ( ERROR " Nothing to build.\n"
    "Please choose the target to build (set the BUILD CMake's"
    " variable to \"MMG\", \"MMG2D\", \"MMGS\" OR"
    " \"MMG3D\")" )

ENDIF()

IF ( BUILD_MMG3D OR BUILD_MMGS )
 SET(BUILD_MMGS3D ON)
ELSE ()
  SET(BUILD_MMGS3D OFF)
ENDIF()

FIND_PACKAGE(Perl)
IF (NOT PERL_FOUND)
  message(WARNING "Perl not found -> we do not generate Fortran Headers")
ENDIF(NOT PERL_FOUND)
############################################################################
#####
#####         Fortran header: libmmgtypesf.h
#####
############################################################################
ADD_EXECUTABLE(genheader ${PROJECT_SOURCE_DIR}/scripts/genheader.c)
TARGET_LINK_LIBRARIES ( genheader PRIVATE ${LIBRARIES}  )
IF (PERL_FOUND)
  GENERATE_FORTRAN_HEADER ( mmg
    ${MMGCOMMON_SOURCE_DIR} libmmgtypes.h
    ${MMGCOMMON_BINARY_DIR}
    ${MMGCOMMON_BINARY_DIR} libmmgtypesf.h
    )
ENDIF (PERL_FOUND)

###############################################################################
#####
#####         libraries
#####
###############################################################################

# Libraries
IF ( NOT (WIN32 OR MINGW) )

  IF(NOT DEFINED M_LIB)
    MESSAGE(STATUS "M_LIB not defined. Searching it")
    FIND_LIBRARY(M_LIB m)
  ENDIF(NOT DEFINED M_LIB)

  MARK_AS_ADVANCED(FORCE M_LIB)
  SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()


############################################################################
#####
#####         Point map
#####
############################################################################

# add point map?
OPTION ( USE_POINTMAP "Use map for point tracking" OFF )

############################################################################
#####
#####         Search for optional libraries: scotch, vtk, linearelasticity
#####
############################################################################
INCLUDE(optional-dependencies)

############################################################################
#####
##### RPATH for MacOSX
#####
############################################################################

SET(CMAKE_MACOSX_RPATH 1)

###############################################################################
#####
#####         Add Mmg options
#####
###############################################################################

INCLUDE(add-mmg-options)

###############################################################################
#####
#####         Build Targets and libmmg if needed
#####
###############################################################################

INCLUDE(build-targets)

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################
CMAKE_DEPENDENT_OPTION (
  ONLY_VERY_SHORT_TESTS "Enable/Disable very short tests" OFF
  "BUILD_TESTING;NOT LONG_TESTS" OFF )
MARK_AS_ADVANCED(ONLY_VERY_SHORT_TESTS)

CMAKE_DEPENDENT_OPTION (
  LONG_TESTS "Enable/Disable time consuming tests for the project" OFF
  "BUILD_TESTING" OFF )
MARK_AS_ADVANCED(LONG_TESTS)

IF( BUILD_TESTING )

  IF ( NOT ONLY_VERY_SHORT_TESTS )
    # Get the continuous integration tests
    INCLUDE(LoadCiTests)
  ENDIF ( )

  # Set the continuous integration options
  SET ( MMG_CI_TESTS ${CI_DIR}/mmg )

  # Copy config file and custom file in the good directories
  CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/cmake/config/CTestCustom.cmake
    ./CTestCustom.cmake)

  ENABLE_TESTING()
  INCLUDE(CTest)

  IF ( ( MMG3D_CI OR MMGS_CI ) AND NOT ONLY_VERY_SHORT_TESTS )
    # Add common tests for mmgs/3d appli
    INCLUDE( mmg_tests )

  ENDIF ( )

ENDIF ( BUILD_TESTING )

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
include(doxygen)

###############################################################################
#####
#####         Create Package
#####
###############################################################################
include(packages)
