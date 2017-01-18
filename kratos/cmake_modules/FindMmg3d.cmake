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
# This package define the MMG3D_INCLUDE_DIR and the MMG3D_LIBRARY variables.  To
# link with the mmg3d library using CMake add the following lines to your
# CMakeLists.txt:
#
# INCLUDE(FindMmg3d.cmake)
#
# INCLUDE_DIRECTORIES(${MMG3D_INCLUDE_DIR})
#
# TARGET_LINK_LIBRARIES( ${YOUR_TARGET} ${MMG3D_LIBRARY})


IF ((NOT WIN32) AND (NOT WIN64))
  SET ( MMG3D_INCLUDE_DIR MMG3D_INCLUDE_DIR-NOTFOUND )
  SET ( MMG3D_LIBRARY MMG3D_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(MMG3D_INCLUDE_DIR
  NAMES mmg/mmg3d/libmmg3d.h
  HINTS ${MMG3D_INCLUDE_DIR}
  $ENV{MMG3D_INCLUDE_DIR}
  $ENV{HOME}/include/
  ${MMG_DIR}/include/
  $ENV{MMG_DIR}/include/
  PATH_SUFFIXES
  DOC "Directory of mmg3d Headers")

# Check for mmg3d library
FIND_LIBRARY(MMG3D_LIBRARY
  NAMES mmg3d mmg3d${MMG3D_LIB_SUFFIX}
  HINTS ${MMG3D_LIBRARY}
  $ENV{MMG3D_LIBRARY}
  $ENV{HOME}/lib
  ${MMG_DIR}/lib
  $ENV{MMG_DIR}/lib
  DOC "The mmg3d library"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MMG3D DEFAULT_MSG
  MMG3D_INCLUDE_DIR MMG3D_LIBRARY)

IF ( NOT MMG3D_FOUND )
   MESSAGE (WARNING "\nMMG3D not found: you can set the 'MMG_DIR' environment"
   " variable or the 'MMG_DIR' CMake variable to your MMG directory path"
   " to help us to fill the MMG3D library paths.\n" )
ENDIF ()
