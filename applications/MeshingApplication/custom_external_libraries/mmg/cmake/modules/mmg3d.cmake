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
## Compilation of mmg3d executable, libraries and tests
##
## =============================================================================

SET(MMG3D_SOURCE_DIR      ${PROJECT_SOURCE_DIR}/src/mmg3d)
SET(MMG3D_BINARY_DIR      ${PROJECT_BINARY_DIR}/src/mmg3d)
SET(MMG3D_SHRT_INCLUDE    mmg/mmg3d )
SET(MMG3D_INCLUDE         ${PROJECT_BINARY_DIR}/include/${MMG3D_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY ${MMG3D_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmg3df.h
#####
############################################################################


if (PERL_FOUND)
  GENERATE_FORTRAN_HEADER ( mmg3d
    ${MMG3D_SOURCE_DIR} libmmg3d.h
    mmg/common
    ${MMG3D_BINARY_DIR} libmmg3df.h
    )
endif (PERL_FOUND)

############################################################################
#####
#####         Choose type of insertion: delaunay kernel or patterns
#####
############################################################################
# Should we use patterns for insertion istead of delaunay kernel
IF ( NOT MMG_PATTERN )
  MESSAGE(STATUS "Vertex insertion by delaunay kernel")
ELSE()
  MESSAGE(STATUS "Vertex insertion by patterns")
  SET(CMAKE_C_FLAGS "-DMMG_PATTERN ${CMAKE_C_FLAGS}")
ENDIF()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Library files
FILE(
  GLOB
  mmg3d_library_files
  ${MMG3D_SOURCE_DIR}/*.c
  ${MMGCOMMON_SOURCE_DIR}/*.c
  ${MMG3D_SOURCE_DIR}/*.h
  ${MMGCOMMON_SOURCE_DIR}/*.h
  ${MMG3D_SOURCE_DIR}/inoutcpp_3d.cpp
  )
LIST(REMOVE_ITEM mmg3d_library_files
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  ${MMGCOMMON_SOURCE_DIR}/apptools.c
  )

IF ( VTK_FOUND AND NOT USE_VTK MATCHES OFF )
  LIST(APPEND  mmg3d_library_files
    ${MMGCOMMON_SOURCE_DIR}/vtkparser.cpp)
ENDIF ( )

FILE(
  GLOB
  mmg3d_main_file
  ${MMG3D_SOURCE_DIR}/mmg3d.c
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
ENDIF()


############################################################################
#####
#####         Compile mmg3d libraries
#####
############################################################################

# Compile static library
IF ( LIBMMG3D_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_a STATIC copy_3d_headers
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# Compile shared library
IF ( LIBMMG3D_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_so SHARED copy_3d_headers
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# mmg3d header files needed for library
#
# Remark: header installation would need to be cleaned, for now, to allow
# independent build of each project and because mmgs and mmg2d have been added
# to mmg3d without rethinking the install architecture, the header files that
# are common between codes are copied in all include directories (mmg/,
# mmg/mmg3d/, mmg/mmgs/, mmg/mmg2d/).  they are also copied in build directory
# to enable library call without installation.
SET( mmg3d_headers
  ${MMG3D_SOURCE_DIR}/mmg3d_export.h
  ${MMG3D_SOURCE_DIR}/libmmg3d.h
  )

IF ( PERL_FOUND )
  LIST ( APPEND mmg3d_headers   ${MMG3D_BINARY_DIR}/libmmg3df.h )
ENDIF()

IF ( MMG_INSTALL_PRIVATE_HEADERS )
  LIST ( APPEND mmg3d_headers
    ${MMG3D_SOURCE_DIR}/libmmg3d_private.h
    ${MMG3D_SOURCE_DIR}/inlined_functions_3d_private.h
    ${MMG3D_SOURCE_DIR}/mmg3dexterns_private.h
    ${MMG3D_SOURCE_DIR}/PRoctree_3d_private.h
    )
ENDIF()

# install man pages
INSTALL(FILES ${PROJECT_SOURCE_DIR}/doc/man/mmg3d.1.gz DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg3d COMPONENT headers)

IF ( MMG_INSTALL_PRIVATE_HEADERS )
  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMG3D_SOURCE_DIR} inlined_functions_3d_private ${MMG3D_INCLUDE} 3d)

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${MMG3D_SOURCE_DIR} PRoctree_3d_private ${MMG3D_INCLUDE} 3d)

  LIST ( APPEND tgt_opt_list copy3d_inlined_functions_3d_private
    copy3d_PRoctree_3d_private )

ENDIF()

# Copy header files in project directory at build step
COPY_HEADERS_AND_CREATE_TARGET ( ${MMG3D_SOURCE_DIR} ${MMG3D_BINARY_DIR} ${MMG3D_INCLUDE} 3d )

###############################################################################
#####
#####         Compile MMG3D executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}3d copy_3d_headers
  "${mmg3d_library_files}" "${mmg3d_main_file}" )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

SET(MMG3D_CI_TESTS ${CI_DIR}/mmg3d )
SET(MMG_CI_TESTS ${CI_DIR}/mmg )

##-------------------------------------------------------------------##
##-------------- Library examples and APIs      ---------------------##
##-------------------------------------------------------------------##
IF ( TEST_LIBMMG3D )
  # Build executables for library examples and add related tests if needed
  INCLUDE(libmmg3d_tests)
ENDIF()

##-------------------------------------------------------------------##
##----------------------- Test Mmg3d executable ---------------------##
##-------------------------------------------------------------------##
IF ( BUILD_TESTING )

  # Add runtime that we want to test for mmg3d
  IF ( MMG3D_CI )

    IF ( LONG_TESTS )
      # Run some tests twice with the output of the previous test as input
      OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )
      MARK_AS_ADVANCED(RUN_AGAIN)
    ELSE ( )
      SET ( RUN_AGAIN OFF )
    ENDIF ( )

    SET ( EXECUT_MMG3D      $<TARGET_FILE:${PROJECT_NAME}3d> )
    SET ( SHRT_EXECUT_MMG3D ${PROJECT_NAME}3d )

    IF ( ONLY_VERY_SHORT_TESTS )
      # Add tests that doesn't require to download meshes
      SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )

      ADD_TEST(NAME mmg3d_very_short COMMAND ${EXECUT_MMG3D}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_a-cube.o"
        )
    ELSE ( )
      # Add mmg3d tests that require to download meshes
      INCLUDE( mmg3d_tests )

    ENDIF ( )

  ENDIF ( MMG3D_CI )

ENDIF ( BUILD_TESTING )
