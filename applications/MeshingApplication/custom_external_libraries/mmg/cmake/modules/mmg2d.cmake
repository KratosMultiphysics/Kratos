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
## Compilation of mmg2d executable, libraries and tests
##
## =============================================================================

SET(MMG2D_SOURCE_DIR      ${PROJECT_SOURCE_DIR}/src/mmg2d)
SET(MMG2D_BINARY_DIR      ${PROJECT_BINARY_DIR}/src/mmg2d)
SET(MMG2D_SHRT_INCLUDE    mmg/mmg2d )
SET(MMG2D_INCLUDE         ${PROJECT_BINARY_DIR}/include/${MMG2D_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY  ${MMG2D_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmg2df.h
#####
############################################################################

if (PERL_FOUND)
  GENERATE_FORTRAN_HEADER ( mmg2d
    ${MMG2D_SOURCE_DIR} libmmg2d.h
    mmg/common
    ${MMG2D_BINARY_DIR} libmmg2df.h
    )
endif(PERL_FOUND)
###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Source files
FILE(
  GLOB
  mmg2d_library_files
  ${MMG2D_SOURCE_DIR}/*.c
  ${MMGCOMMON_SOURCE_DIR}/*.c
  ${MMG2D_SOURCE_DIR}/*.h
  ${MMGCOMMON_SOURCE_DIR}/*.h
  ${MMG2D_SOURCE_DIR}/inoutcpp_2d.cpp
  )

LIST(REMOVE_ITEM mmg2d_library_files
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  ${MMGCOMMON_SOURCE_DIR}/apptools.c
  ${REMOVE_FILE} )

IF ( VTK_FOUND AND NOT USE_VTK MATCHES OFF )
  LIST(APPEND  mmg2d_library_files
    ${MMGCOMMON_SOURCE_DIR}/vtkparser.cpp )
ENDIF ( )

FILE(
  GLOB
  mmg2d_main_file
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  ${MMGCOMMON_SOURCE_DIR}/apptools.c
  )

############################################################################
#####
#####         Elastic
#####
############################################################################

IF( ELAS_FOUND AND NOT USE_ELAS MATCHES OFF )
  # Set flags for building test program
  INCLUDE_DIRECTORIES(AFTER ${ELAS_INCLUDE_DIR})

  SET(CMAKE_REQUIRED_INCLUDES ${ELAS_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${ELAS_LIBRARY})

  SET(CMAKE_C_FLAGS "-DUSE_ELAS ${CMAKE_C_FLAGS}")
  MESSAGE(STATUS
    "Compilation with the Elas library: ${ELAS_LIBRARY} ")
  SET( LIBRARIES ${ELAS_LINK_FLAGS} ${LIBRARIES})
  SET( LIBRARIES ${ELAS_LIBRARY} ${LIBRARIES})

ENDIF ( )

############################################################################
#####
#####         Compile mmg2d libraries
#####
############################################################################
# mmg2d header files needed for library
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

IF ( PERL_FOUND )
  LIST ( APPEND mmg2d_headers   ${MMG2D_BINARY_DIR}/libmmg2df.h )
ENDIF()

IF ( MMG_INSTALL_PRIVATE_HEADERS )
  LIST ( APPEND mmg2d_headers
    ${MMG2D_SOURCE_DIR}/libmmg2d_private.h
    ${MMG2D_SOURCE_DIR}/mmg2dexterns_private.h
    )
ENDIF()

# install man pages
INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmg2d.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmg2d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg2d COMPONENT headers )

# Copy header files in project directory at build step
COPY_HEADERS_AND_CREATE_TARGET ( ${MMG2D_SOURCE_DIR} ${MMG2D_BINARY_DIR} ${MMG2D_INCLUDE} 2d )

# Compile static library
IF ( LIBMMG2D_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}2d_a STATIC copy_2d_headers
    "${mmg2d_library_files}" ${PROJECT_NAME}2d )
ENDIF()

# Compile shared library
IF ( LIBMMG2D_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}2d_so SHARED copy_2d_headers
    "${mmg2d_library_files}" ${PROJECT_NAME}2d )
ENDIF()

###############################################################################
#####
#####         Compile MMG2D executable
#####
###############################################################################

ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}2d copy_2d_headers
  "${mmg2d_library_files}" "${mmg2d_main_file}" )


###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

SET(MMG2D_CI_TESTS ${CI_DIR}/mmg2d )

##-------------------------------------------------------------------##
##-------------- Library examples and APIs      ---------------------##
##-------------------------------------------------------------------##
IF ( TEST_LIBMMG2D )
  # Build library executables and add library tests if needed
  INCLUDE(libmmg2d_tests)
ENDIF ( )

##-------------------------------------------------------------------##
##----------------------- Test Mmg2d executable ---------------------##
##-------------------------------------------------------------------##
IF ( BUILD_TESTING )

  # Add runtime that we want to test for mmg2d
  IF ( MMG2D_CI )

    SET ( EXECUT_MMG2D      $<TARGET_FILE:${PROJECT_NAME}2d> )

    IF ( ONLY_VERY_SHORT_TESTS )
      # Add tests that doesn't require to download meshes
      SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )

      ADD_TEST(NAME mmg2d_very_short COMMAND ${EXECUT_MMG2D}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/init.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_a-init.o"
        )
    ELSE ( )
      # Add mmg2d tests that require to download meshes
      INCLUDE( mmg2d_tests )

    ENDIF ( )

  ENDIF( MMG2D_CI )

ENDIF ( BUILD_TESTING )
