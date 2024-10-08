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

# Build Doxygen doc if needed
OPTION ( BUILD_DOC "Enable/Disable Doxygen documentation target" OFF )
IF ( BUILD_DOC )
  FIND_PACKAGE(Doxygen)
  IF(DOXYGEN_FOUND)

    # MMG Documentation
    SET( MMG_DOC_DIR ${PROJECT_BINARY_DIR}/doc/ CACHE PATH
      "Path toward generated Doxygen doc.")
    mark_as_advanced(MMG_DOC_DIR)

    FILE(MAKE_DIRECTORY ${MMG_DOC_DIR})

    CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    ADD_CUSTOM_TARGET(doc
      COMMAND ${DOXYGEN_EXECUTABLE}
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${MMG_DOC_DIR}
      COMMENT "Generating Mmg documentation with Doxygen. Open up the
    ${MMG_DOC_DIR}/index.html file to see
    it" VERBATIM )

  ENDIF ( DOXYGEN_FOUND )
ENDIF ( )
