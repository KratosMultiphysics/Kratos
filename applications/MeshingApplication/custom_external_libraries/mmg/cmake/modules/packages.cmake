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

###############################################################################
#####
#####         Create Package
#####
###############################################################################

INCLUDE(InstallRequiredSystemLibraries)
SET( CPACK_ARCHIVE_COMPONENT_INSTALL ON )
SET( CPACK_COMPONENT_ALL appli )
SET( CPACK_SOURCE_GENERATOR TGZ )
SET( CPACK_GENERATOR TGZ )

SET(CPACK_PACKAGE_NAME ${PROJECT_NAME})
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG: 2d, surface and 3d remeshers")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey, Charles Dapogny,"
  " Algiane Froehly")
SET(CPACK_PACKAGE_CONTACT "contact@mmgtools.org")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION ${CMAKE_RELEASE_VERSION})
SET(CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_RELEASE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${CMAKE_RELEASE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${CMAKE_RELEASE_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_INSTALL_DIRECTORY}-${CMAKE_SYSTEM})
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}")
SET(CPACK_OUTPUT_FILE_PREFIX "../archives")

INCLUDE(CPack)
