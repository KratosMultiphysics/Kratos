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
  SET ( SCOTCH_INCLUDE_DIR SCOTCH_INCLUDE_DIR-NOTFOUND )
  SET ( SCOTCH_LIBRARY SCOTCH_LIBRARY-NOTFOUND )
  SET ( SCOTCHERR_LIBRARY SCOTCHERR_LIBRARY-NOTFOUND )
ENDIF()

FIND_PATH(SCOTCH_INCLUDE_DIR
  NAMES scotch.h
  HINTS ${SCOTCH_INCLUDE_DIR}
  $ENV{SCOTCH_INCLUDE_DIR}
  ${SCOTCH_DIR}/include
  $ENV{SCOTCH_DIR}/include
  PATH_SUFFIXES scotch
  DOC "Directory of SCOTCH Header")

# Check for scotch
FIND_LIBRARY(SCOTCH_LIBRARY
  NAMES scotch scotch${SCOTCH_LIB_SUFFIX}
  HINTS ${SCOTCH_LIBRARY}
  $ENV{SCOTCH_LIBRARY}
  ${SCOTCH_DIR}/lib
  $ENV{SCOTCH_DIR}/lib
  DOC "The SCOTCH library"
  )

FIND_LIBRARY(SCOTCHERR_LIBRARY
  NAMES scotcherr scotcherr${SCOTCH_LIB_SUFFIX}
  HINTS ${SCOTCHERR_LIBRARY}
  $ENV{SCOTCHERR_LIBRARY}
  ${SCOTCH_DIR}/lib
  $ENV{SCOTCH_DIR}/lib
  DOC "The SCOTCH-ERROR library"
  )

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(SCOTCH DEFAULT_MSG
				  SCOTCH_INCLUDE_DIR SCOTCH_LIBRARY SCOTCHERR_LIBRARY)
IF ((NOT WIN32) AND (NOT WIN64))
  MARK_AS_ADVANCED(SCOTCH_INCLUDE_DIR SCOTCH_LIBRARY SCOTCHERR_LIBRARY)
ENDIF()
