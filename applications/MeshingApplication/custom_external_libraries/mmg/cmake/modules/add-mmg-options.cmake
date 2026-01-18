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

#------------------------------- static or shared libs
FUNCTION (INVERT_BOOL OUTVAR INVAR)
  IF(${INVAR})
    SET(${OUTVAR} OFF PARENT_SCOPE)
  ELSE()
    SET(${OUTVAR} ON PARENT_SCOPE)
  ENDIF()
ENDFUNCTION()

OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)
INVERT_BOOL("BUILD_STATIC_LIBS" ${BUILD_SHARED_LIBS})
IF ( (${BUILD_STATIC_LIBS} EQUAL ON) AND NOT CMAKE_POSITION_INDEPENDENT_CODE)
  SET(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

#------------------------------- mmg2d
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_STATIC "Compile the mmg2d static library" ON
  "BUILD_MMG2D;BUILD_STATIC_LIBS" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_SHARED "Compile the mmg2d dynamic library" ON
  "BUILD_MMG2D;BUILD_SHARED_LIBS" OFF )
IF (LIBMMG2D_STATIC OR LIBMMG2D_SHARED)
  SET(LIBMMG2D_INTERNAL ON )
ELSE()
  SET(LIBMMG2D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG2D  "Build mmg2d library examples and tests" OFF
  "BUILD_MMG2D;LIBMMG2D_INTERNAL" OFF)

CMAKE_DEPENDENT_OPTION(
  MMG2D_CI "Enable/Disabletest_execution for mmg2d" ON
  "BUILD_TESTING;BUILD_MMG2D" OFF)

#------------------------------- mmgs
CMAKE_DEPENDENT_OPTION ( LIBMMGS_STATIC "Compile the mmgs static library" ON
  "BUILD_MMGS;BUILD_STATIC_LIBS" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMGS_SHARED "Compile the mmgs dynamic library" ON
  "BUILD_MMGS;BUILD_SHARED_LIBS" OFF)
IF (LIBMMGS_STATIC OR LIBMMGS_SHARED)
  SET(LIBMMGS_INTERNAL ON )
ELSE()
  SET(LIBMMGS_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMGS "Build mmgs library examples and tests" OFF
  "LIBMMGS_INTERNAL;BUILD_MMGS" OFF
  )

CMAKE_DEPENDENT_OPTION(
  MMGS_CI "Enable/Disable test_execution for mmgs" ON
  "BUILD_TESTING;BUILD_MMGS" OFF )

#------------------------------- mmg3d
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_STATIC "Compile the mmg3d static library" ON
  "BUILD_MMG3D;BUILD_STATIC_LIBS" OFF )
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_SHARED "Compile the mmg3d dynamic library" ON
  "BUILD_MMG3D;BUILD_SHARED_LIBS" OFF )
IF (LIBMMG3D_STATIC OR LIBMMG3D_SHARED)
  SET(LIBMMG3D_INTERNAL ON )
ELSE()
  SET(LIBMMG3D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG3D "Build mmg3d library examples and tests" OFF
  "BUILD_MMG3D;LIBMMG3D_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG3D_CI "Enable/Disable test execution for mmg3d" ON
  "BUILD_TESTING;BUILD_MMG3D" OFF)

CMAKE_DEPENDENT_OPTION ( MMG_PATTERN
  "If ON, insertion by patterns, otherwise insertion by delaunay kernel"
  OFF
  "BUILD_MMG3D" OFF)
#-- Remove the next line to have this option visible in basic cmake mode
MARK_AS_ADVANCED(MMG_PATTERN)

#---------------------------- mmg library gathering mmg2d, mmgs and mmg3d libs
CMAKE_DEPENDENT_OPTION ( LIBMMG_STATIC
  "Compile the mmg static library (mmg2d + mmgs + mmg3d)" ON
  "BUILD_MMG;BUILD_STATIC_LIBS" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG_SHARED
  "Compile the mmg dynamic library (mmg2d + mmgs + mmg3d)" ON
  "BUILD_MMG;BUILD_SHARED_LIBS" OFF )
IF (LIBMMG_STATIC OR LIBMMG_SHARED)
  SET(LIBMMG_INTERNAL ON )
ELSE()
  SET(LIBMMG_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG "Build mmg library examples and tests" OFF
  "BUILD_MMG;LIBMMG_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG_CI "Enable/Disable test execution for mmg library" ON
  "BUILD_TESTING;BUILD_MMG" OFF)

#---------------------- install provate headers for library built on top of Mmg
OPTION ( MMG_INSTALL_PRIVATE_HEADERS
  "Install private headers of Mmg (use at your own risk)" OFF )
