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


GENERATE_FORTRAN_HEADER ( mmg3d
  ${MMG3D_SOURCE_DIR} libmmg3d.h
  ${MMG3D_SHRT_INCLUDE}
  ${MMG3D_BINARY_DIR} libmmg3df.h
  )

############################################################################
#####
#####         Choose type of insertion: delaunay kernel or patterns
#####
############################################################################
# Should we use patterns for insertion istead of delaunay kernel
IF ( NOT PATTERN )
  MESSAGE(STATUS "Vertex insertion by delaunay kernel")
ELSE()
  MESSAGE(STATUS "Vertex insertion by patterns")
  SET(CMAKE_C_FLAGS "-DPATTERN ${CMAKE_C_FLAGS}")
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
  ${COMMON_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM mmg3d_library_files
  ${MMG3D_SOURCE_DIR}/${PROJECT_NAME}3d.c
)
FILE(
  GLOB
  mmg3d_main_file
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  )

############################################################################
#####
#####         Elastic
#####
############################################################################

IF( USE_ELAS )
  # Set flags for building test program
  INCLUDE_DIRECTORIES(${ELAS_INCLUDE_DIR})

  SET(CMAKE_REQUIRED_INCLUDES ${ELAS_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${ELAS_LIBRARY})

  SET(CMAKE_C_FLAGS "-DUSE_ELAS ${CMAKE_C_FLAGS}")
  MESSAGE(STATUS
    "Compilation with the Elas library: ${ELAS_LIBRARY} ")
  SET( LIBRARIES ${ELAS_LINK_FLAGS} ${LIBRARIES})
  SET( LIBRARIES ${ELAS_LIBRARY} ${LIBRARIES})
ENDIF()

IF (ELAS_NOTFOUND)
  MESSAGE ( WARNING "Elas is a library to solve the linear elasticity "
    "problem (see https://github.com/ISCDtoolbox/LinearElasticity to"
    " download it). "
    "This library is needed to use the lagrangian motion option. "
    "If you have already installed Elas and want to use it, "
    "please set the CMake variable or environment variable ELAS_DIR "
    "to your Elas directory.")
ENDIF ( )

############################################################################
#####
#####         Compile mmg3d libraries
#####
############################################################################

# Compile static library
IF ( LIBMMG3D_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_a STATIC
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# Compile shared library
IF ( LIBMMG3D_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}3d_so SHARED
    "${mmg3d_library_files}" ${PROJECT_NAME}3d )
ENDIF()

# mmg3d header files needed for library
SET( mmg3d_headers
  ${MMG3D_SOURCE_DIR}/libmmg3d.h
  ${MMG3D_BINARY_DIR}/libmmg3df.h
  ${COMMON_SOURCE_DIR}/libmmgtypes.h
  ${COMMON_BINARY_DIR}/libmmgtypesf.h
  )

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg3d COMPONENT headers)

COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${MMG3D_BINARY_DIR} ${MMG3D_INCLUDE} 3d )

# Copy header files in project directory at configuration step
# (generated file don't exists yet or are outdated)
FILE(INSTALL  ${mmg3d_headers} DESTINATION ${MMG3D_INCLUDE}
  PATTERN "libmmg*f.h"  EXCLUDE)

############################################################################
#####
#####         Compile program to test library
#####
############################################################################
SET(MMG3D_CI_TESTS ${CI_DIR}/mmg3d )
SET(MMG_CI_TESTS ${CI_DIR}/mmg )

IF ( TEST_LIBMMG3D )
  INCLUDE(cmake/testing/libmmg3d_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile MMG3D executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}3d
  "${mmg3d_library_files}" ${mmg3d_main_file} )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg3d
  IF ( MMG3D_CI )

    SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )
    FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )

    IF ( LONG_TESTS )
      # Run some tests twice with the output of the previous test as input
      OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )
      MARK_AS_ADVANCED(RUN_AGAIN)
    ELSE ( )
      SET ( RUN_AGAIN OFF )
    ENDIF ( )

    ADD_EXEC_TO_CI_TESTS ( ${PROJECT_NAME}3d EXECUT_MMG3D )
    SET ( LISTEXEC_MMG ${EXECUT_MMG3D} )

    IF ( TEST_LIBMMG3D )
      SET(LIBMMG3D_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_a)
      SET(LIBMMG3D_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_b)
      SET(LIBMMG3D_EXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example1)
      SET(LIBMMG3D_EXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example2)
      SET(LIBMMG3D_EXEC4   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example4)
      SET(LIBMMG3D_EXEC5   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example5)
      SET(LIBMMG3D_EXEC6   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example6_io)
      SET(TEST_API3D_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api3d_0)

      ADD_TEST(NAME libmmg3d_example0_a COMMAND ${LIBMMG3D_EXEC0_a}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_a-cube.o"
        )
      ADD_TEST(NAME libmmg3d_example0_b COMMAND ${LIBMMG3D_EXEC0_b}
       "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_b.o.mesh"
        )
      ADD_TEST(NAME libmmg3d_example1   COMMAND ${LIBMMG3D_EXEC1}
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1.o.mesh"
        )
      ADD_TEST(NAME libmmg3d_example2   COMMAND ${LIBMMG3D_EXEC2}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example2/2spheres.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1-2spheres_1.o"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_1-2spheres_2.o"
        )
      IF ( USE_ELAS )
        ADD_TEST(NAME libmmg3d_example4   COMMAND ${LIBMMG3D_EXEC4}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/LagrangianMotion_example0/tinyBoxt"
          "${CTEST_OUTPUT_DIR}/libmmg3d_LagrangianMotion_0-tinyBoxt.o"
          )
      ENDIF ()
      ADD_TEST(NAME libmmg3d_example5   COMMAND ${LIBMMG3D_EXEC5}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/IsosurfDiscretization_example0/test"
        "${CTEST_OUTPUT_DIR}/libmmg3d-IsosurfDiscretization_0-test.o"
        )
      ADD_TEST(NAME libmmg3d_example6_io_0   COMMAND ${LIBMMG3D_EXEC6}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_io_6-naca.o" "0"
       )
      ADD_TEST(NAME libmmg3d_example6_io_1   COMMAND ${LIBMMG3D_EXEC6}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_io_6-naca.o" "1"
       )
      ADD_TEST(NAME test_api3d_0   COMMAND ${TEST_API3D_EXEC0}
        "${MMG3D_CI_TESTS}/API_tests/2dom.mesh"
        "${CTEST_OUTPUT_DIR}/test_API3d.o"
       )

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG3D_EXECFORTRAN_a  ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_a)
        SET(LIBMMG3D_EXECFORTRAN_b  ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_b)
        SET(LIBMMG3D_EXECFORTRAN_IO ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_io)
        SET(TEST_API3D_FORTRAN_EXEC0 ${EXECUTABLE_OUTPUT_PATH}/test_api3d_fortran_0)


        ADD_TEST(NAME libmmg3d_fortran_a  COMMAND ${LIBMMG3D_EXECFORTRAN_a}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0_fortran/example0_a/cube.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d-Adaptation_Fortran_0_a-cube.o"
          )
        ADD_TEST(NAME libmmg3d_fortran_b  COMMAND ${LIBMMG3D_EXECFORTRAN_b}
          "${CTEST_OUTPUT_DIR}/libmmg3d-Adaptation_Fortran_0_b-cube.o"
          )
        ADD_TEST(NAME libmmg3d_fortran_io_0   COMMAND ${LIBMMG3D_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d_Fortran_io-torus.o" "0"
          )
        ADD_TEST(NAME libmmg3d_fortran_io_1   COMMAND ${LIBMMG3D_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/io_multisols_example6/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmg3d_Fortran_io-torus.o" "1"
         )
       ADD_TEST(NAME test_api3d_fortran_0   COMMAND ${TEST_API3D_FORTRAN_EXEC0}
         "${MMG3D_CI_TESTS}/API_tests/2dom.mesh"
         "${CTEST_OUTPUT_DIR}/test_API3d.o"
         )

      ENDIF()

    ENDIF ( TEST_LIBMMG3D )

    IF ( ONLY_VERY_SHORT_TESTS )
      SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )

      ADD_TEST(NAME mmg3d_very_short COMMAND ${EXECUT_MMG3D}
        "${PROJECT_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmg3d_Adaptation_0_a-cube.o"
        )
    ELSE ( )

      # Add more tests
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmg3d_tests.cmake )
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmg_tests.cmake )
    ENDIF ( )

  ENDIF ( MMG3D_CI )

ENDIF ( BUILD_TESTING )
