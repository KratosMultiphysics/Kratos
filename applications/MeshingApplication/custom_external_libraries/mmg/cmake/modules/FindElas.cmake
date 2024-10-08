## =============================================================================
##  This file is part of the mmg software package for the tetrahedral
##  mesh modification.
##**  Copyright (c) Bx INP/Inria/UBordeaux/UPMC, 2004- .
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

IF ((NOT WIN32) AND (NOT WIN64))
  SET ( ELAS_INCLUDE_DIR ELAS_INCLUDE_DIR-NOTFOUND )
  SET ( ELAS_LIBRARY ELAS_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(ELAS_INCLUDE_DIR
  NAMES elastic.h
  HINTS ${ELAS_INCLUDE_DIR}
  $ENV{ELAS_INCLUDE_DIR}
  $ENV{HOME}/include
  ${ELAS_DIR}/include
  $ENV{ELAS_DIR}/include
  ${ELAS_DIR}/sources
  $ENV{ELAS_DIR}/sources
  PATH_SUFFIXES ELAS
  DOC "Directory of ELAS Header")

# Check for elas library (search in common locations in order to prevent
# modifications in ELAS installation)
FIND_LIBRARY(ELAS_LIBRARY
  NAMES Elas Elas${ELAS_LIB_SUFFIX}
  HINTS ${ELAS_LIBRARY}
  $ENV{ELAS_LIBRARY}
  $ENV{HOME}/lib
  ${ELAS_DIR}/lib
  $ENV{ELAS_DIR}/lib
  DOC "The ELAS library"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ELAS DEFAULT_MSG
  ELAS_INCLUDE_DIR ELAS_LIBRARY)
IF ((NOT WIN32) AND (NOT WIN64))
  MARK_AS_ADVANCED(ELAS_INCLUDE_DIR ELAS_LIBRARY)
ENDIF()
