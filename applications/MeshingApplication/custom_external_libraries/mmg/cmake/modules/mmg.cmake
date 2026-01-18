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

## =============================================================================
##
## Compilation of mmg library and tests
##
## =============================================================================

# Note: mmg3d.cmake, mmg2d.cmake and mmgs.cmake must have been called before
# (thus the MMG*_SOURCE_DIR variables are already setted)

############################################################################
#####
#####         Compile mmg libraries
#####
############################################################################

# Compile static library
FILE(
  GLOB
  mmg_library_files
  ${MMG2D_SOURCE_DIR}/inoutcpp_2d.cpp
  ${MMG3D_SOURCE_DIR}/inoutcpp_3d.cpp
  ${MMGS_SOURCE_DIR}/inoutcpp_s.cpp
  ${MMG2D_SOURCE_DIR}/*.c
  ${MMG3D_SOURCE_DIR}/*.c
  ${MMGS_SOURCE_DIR}/*.c
  ${MMGCOMMON_SOURCE_DIR}/*.c
  ${MMG2D_SOURCE_DIR}/*.h
  ${MMG3D_SOURCE_DIR}/*.h
  ${MMGS_SOURCE_DIR}/*.h
  ${MMGCOMMON_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM mmg_library_files
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  ${MMGS_SOURCE_DIR}/mmgs.c
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  ${REMOVE_FILE} )

IF ( VTK_FOUND AND NOT USE_VTK MATCHES OFF)
  LIST(APPEND  mmg_library_files
    ${MMGCOMMON_SOURCE_DIR}/vtkparser.cpp )
ENDIF ( )

IF ( LIBMMG_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_a  STATIC
    copy_mmg_headers
    "${mmg_library_files}"
    ${PROJECT_NAME} )
ENDIF()

# Compile shared library
IF ( LIBMMG_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_so  SHARED
    copy_mmg_headers
    "${mmg_library_files}"
    ${PROJECT_NAME} )
ENDIF()

IF ( LIBMMG_STATIC OR LIBMMG_SHARED )
  # mmg header files needed for library
  #
  # Remark: header installation would need to be cleaned, for now, to allow
  # independent build of each project and because mmgs and mmg2d have been added
  # to mmg3d without rethinking the install architecture, the header files that
  # are common between codes are copied in all include directories (mmg/,
  # mmg/mmg3d/, mmg/mmgs/, mmg/mmg2d/).  they are also copied in build directory
  # to enable library call without installation.
  SET( mmg2d_headers
    ${MMG2D_SOURCE_DIR}/mmg2d_export.h
    ${MMG2D_SOURCE_DIR}/libmmg2d.h
    )
  SET( mmg3d_headers
    ${MMG3D_SOURCE_DIR}/mmg3d_export.h
    ${MMG3D_SOURCE_DIR}/libmmg3d.h
    )
  SET( mmgs_headers
    ${MMGS_SOURCE_DIR}/mmgs_export.h
    ${MMGS_SOURCE_DIR}/libmmgs.h
    )
  SET( mmg_headers
     # ${PROJECT_SOURCE_DIR}/src/common/mmg_core_export.h
     ${PROJECT_SOURCE_DIR}/src/mmg/libmmg.h
    )

  IF ( PERL_FOUND )
    LIST ( APPEND mmg3d_headers   ${MMG3D_BINARY_DIR}/libmmg3df.h )
    LIST ( APPEND mmg2d_headers   ${MMG2D_BINARY_DIR}/libmmg2df.h )
    LIST ( APPEND mmgs_headers    ${MMGS_BINARY_DIR}/libmmgsf.h )
    LIST ( APPEND mmg_headers     ${PROJECT_SOURCE_DIR}/src/mmg/libmmgf.h )
  ENDIF()


  SET(MMG2D_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmg2d )
  SET(MMGS_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmgs )
  SET(MMG3D_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmg3d )
  SET(MMG_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg )

  # install man pages
  INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmg2d.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
  INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmg3d.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
  INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmgs.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg)

  IF ( PERL_FOUND )
    FILE(INSTALL ${PROJECT_SOURCE_DIR}/src/mmg/libmmgf.h DESTINATION  ${PROJECT_BINARY_DIR}/include/mmg/)
  ENDIF()

  # Install header files in project directory
  FILE(INSTALL  ${mmg2d_headers}
    DESTINATION ${PROJECT_BINARY_DIR}/include/mmg/mmg2d/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmgs_headers}
    DESTINATION ${PROJECT_BINARY_DIR}/include/mmg/mmgs/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmg3d_headers}
    DESTINATION ${PROJECT_BINARY_DIR}/include/mmg/mmg3d/
    PATTERN "libmmg*f.h"  EXCLUDE)
  FILE(INSTALL  ${mmg_headers}
    DESTINATION ${PROJECT_BINARY_DIR}/include/mmg/
    PATTERN "libmmg*f.h"  EXCLUDE)

  set ( mmg_file_list ${PROJECT_BINARY_DIR}/include/mmg/libmmg.h )
  IF ( PERL_FOUND )
    list ( APPEND mmg_file_list ${PROJECT_BINARY_DIR}/include/mmg/libmmgf.h)
  ENDIF( )

  ADD_CUSTOM_TARGET(copy_mmg_headers ALL
    DEPENDS
    copy_2d_headers copy_s_headers copy_3d_headers ${mmg_file_list}
    )

ENDIF()

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

##-------------------------------------------------------------------##
##-------------- Library examples and APIs      ---------------------##
##-------------------------------------------------------------------##
IF ( TEST_LIBMMG )
  INCLUDE(libmmg_tests)
ENDIF()
