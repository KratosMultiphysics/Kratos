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

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# CTEST_OUTPUT_DIR MUST BE SETTED HERE TO AVOID UNINITIALIZATION FOR LIBRARY
# TESTS DEFS
OPTION ( BUILD_TESTING "Enable/Disable continuous integration" OFF )

IF( BUILD_TESTING )
  SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS
    CACHE PATH "Path toward the tests outputs" )
  MARK_AS_ADVANCED ( CTEST_OUTPUT_DIR )
  FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )
ENDIF ( )

INCLUDE(mmgcommon)

IF ( BUILD_MMG2D )
  INCLUDE(mmg2d)
ENDIF ( )

IF ( BUILD_MMGS )
  INCLUDE(mmgs)
ENDIF ( )

IF ( BUILD_MMG3D )
  INCLUDE(mmg3d)
ENDIF ( )

IF ( BUILD_MMG )
  INCLUDE(mmg)
ENDIF ( )

IF ( MmgTargetsExported )
  set( MMG_CMAKE_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/mmg )

  configure_package_config_file(cmake/config/mmgConfig.cmake.in
    ${PROJECT_BINARY_DIR}/mmgConfig.cmake
    INSTALL_DESTINATION ${MMG_CMAKE_INSTALL_DIR} )

  install(EXPORT MmgTargets
    FILE MmgTargets.cmake
    NAMESPACE Mmg::
    DESTINATION ${MMG_CMAKE_INSTALL_DIR}
    )
  install(FILES ${PROJECT_BINARY_DIR}/mmgConfig.cmake
    DESTINATION ${MMG_CMAKE_INSTALL_DIR} )

  install(FILES
    ${PROJECT_SOURCE_DIR}/cmake/modules/FindSCOTCH.cmake
    ${PROJECT_SOURCE_DIR}/cmake/modules/FindElas.cmake
    DESTINATION ${MMG_CMAKE_INSTALL_DIR} )
ENDIF()
