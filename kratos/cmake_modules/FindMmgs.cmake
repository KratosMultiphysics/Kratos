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
# This package define the MMGS_INCLUDE_DIR and the MMGS_LIBRARY variables.  To
# link with the mmgs library using CMake add the following lines to your
# CMakeLists.txt:
#
# INCLUDE(FindMmgs.cmake)
#
# INCLUDE_DIRECTORIES(${MMGS_INCLUDE_DIR})
#
# TARGET_LINK_LIBRARIES( ${YOUR_TARGET} ${MMGS_LIBRARY})


IF ((NOT WIN32) AND (NOT WIN64))
  SET ( MMGS_INCLUDE_DIR MMGS_INCLUDE_DIR-NOTFOUND )
  SET ( MMGS_LIBRARY MMGS_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(MMGS_INCLUDE_DIR
  NAMES mmg/mmgs/libmmgs.h
  HINTS ${MMGS_INCLUDE_DIR}
  $ENV{MMGS_INCLUDE_DIR}
  $ENV{HOME}/include/
  ${MMG_DIR}/include/
  $ENV{MMG_DIR}/include/
  PATH_SUFFIXES
  DOC "Directory of mmgs Headers")

# Check for mmgs library
FIND_LIBRARY(MMGS_LIBRARY
  NAMES mmgs mmgs${MMGS_LIB_SUFFIX}
  HINTS ${MMGS_LIBRARY}
  $ENV{MMGS_LIBRARY}
  $ENV{HOME}/lib
  ${MMG_DIR}/lib
  $ENV{MMG_DIR}/lib
  DOC "The mmgs library"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(MMGS DEFAULT_MSG
  MMGS_INCLUDE_DIR MMGS_LIBRARY)

IF ( NOT MMGS_FOUND )
   MESSAGE (WARNING "\nMMGS not found: you can set the 'MMG_DIR' environment"
   " variable or the 'MMG_DIR' CMake variable to your MMG directory path"
   " to help us to fill the MMGS library paths.\n" )
ENDIF ()
