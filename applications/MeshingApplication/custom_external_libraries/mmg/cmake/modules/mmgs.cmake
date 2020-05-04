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
## Compilation of mmgs executable, libraries and tests
##
## =============================================================================

SET(MMGS_SOURCE_DIR      ${PROJECT_SOURCE_DIR}/src/mmgs)
SET(MMGS_BINARY_DIR      ${PROJECT_BINARY_DIR}/src/mmgs)
SET(MMGS_SHRT_INCLUDE    mmg/mmgs )
SET(MMGS_INCLUDE         ${PROJECT_BINARY_DIR}/include/${MMGS_SHRT_INCLUDE} )

FILE(MAKE_DIRECTORY ${MMGS_BINARY_DIR})

############################################################################
#####
#####         Fortran header: libmmgsf.h
#####
############################################################################

GENERATE_FORTRAN_HEADER ( mmgs
  ${MMGS_SOURCE_DIR} libmmgs.h
  ${MMGS_SHRT_INCLUDE}
  ${MMGS_BINARY_DIR} libmmgsf.h
  )

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################

# Source files
FILE(
  GLOB
  mmgs_library_files
  ${MMGS_SOURCE_DIR}/*.c
  ${COMMON_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM mmgs_library_files
  ${MMGS_SOURCE_DIR}/mmgs.c
  ${REMOVE_FILE} )
FILE(
  GLOB
  mmgs_main_file
  ${MMGS_SOURCE_DIR}/mmgs.c
  )

############################################################################
#####
#####         Compile mmgs libraries
#####
############################################################################

# Compile static library
IF ( LIBMMGS_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}s_a STATIC
    "${mmgs_library_files}" ${PROJECT_NAME}s )
ENDIF()

# Compile shared library
IF ( LIBMMGS_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}s_so SHARED
    "${mmgs_library_files}" ${PROJECT_NAME}s )
ENDIF()

# mmgs header files needed for library
SET( mmgs_headers
  ${MMGS_SOURCE_DIR}/libmmgs.h
  ${MMGS_BINARY_DIR}/libmmgsf.h
  ${COMMON_SOURCE_DIR}/libmmgtypes.h
  ${COMMON_BINARY_DIR}/libmmgtypesf.h
  )

# Install header files in /usr/local or equivalent
INSTALL(FILES ${mmgs_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmgs COMPONENT headers)

COPY_FORTRAN_HEADER_AND_CREATE_TARGET ( ${MMGS_BINARY_DIR} ${MMGS_INCLUDE} s )

# Copy header files in project directory at configuration step
# (generated file don't exists yet or are outdated)
FILE(INSTALL  ${mmgs_headers} DESTINATION ${MMGS_INCLUDE}
  PATTERN "libmmg*f.h"  EXCLUDE)

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMGS )
  INCLUDE(cmake/testing/libmmgs_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile MMGS executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME}s
  "${mmgs_library_files}" ${mmgs_main_file} )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET(MMGS_CI_TESTS ${CI_DIR}/mmgs )
  SET(MMG_CI_TESTS ${CI_DIR}/mmg )

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmgs
  IF( MMGS_CI )

    SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )
    FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )


    ADD_EXEC_TO_CI_TESTS ( ${PROJECT_NAME}s EXECUT_MMGS )
    SET ( LISTEXEC_MMG ${EXECUT_MMGS} )

    IF ( TEST_LIBMMGS )
      SET(LIBMMGS_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_a)
      SET(LIBMMGS_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_b)
      SET(LIBMMGS_EXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example1)
      SET(LIBMMGS_EXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example2)
      SET(LIBMMGS_EXEC3   ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example3)


      ADD_TEST(NAME libmmgs_example0_a   COMMAND ${LIBMMGS_EXEC0_a}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_0_a-cube.o"
        )
      ADD_TEST(NAME libmmgs_example0_b  COMMAND ${LIBMMGS_EXEC0_b}
        "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_0_b.o.mesh"
        )
      ADD_TEST(NAME libmmgs_example1   COMMAND ${LIBMMGS_EXEC1}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/2spheres"
        "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_1-2spheres_1.o"
        "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_1-2spheres_2.o"
        )
      ADD_TEST(NAME libmmgs_example2   COMMAND ${LIBMMGS_EXEC2}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/teapot"
        "${CTEST_OUTPUT_DIR}/libmmgs-IsosurfDiscretization_0-teapot.o"
        )
      ADD_TEST(NAME libmmgs_example3_io_0   COMMAND ${LIBMMGS_EXEC3}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmgs_io_3-naca.o" "0"
       )
      ADD_TEST(NAME libmmgs_example3_io_1   COMMAND ${LIBMMGS_EXEC3}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
        "${CTEST_OUTPUT_DIR}/libmmgs_io_3-naca.o" "1"
       )


      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMGS_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_a)
        SET(LIBMMGS_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_b)
        SET(LIBMMGS_EXECFORTRAN_IO ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_io)

        ADD_TEST(NAME libmmgs_fortran_a   COMMAND ${LIBMMGS_EXECFORTRAN_a}
          "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_a/cube.mesh"
          "${CTEST_OUTPUT_DIR}/libmmgs-Adaptation_Fortran_0_a-cube.o"
          )
        ADD_TEST(NAME libmmgs_fortran_b   COMMAND ${LIBMMGS_EXECFORTRAN_b}
          "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_Fortran_0_b.o"
          )
        ADD_TEST(NAME libmmgs_fortran_io_0   COMMAND ${LIBMMGS_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmgs_Fortran_io-torus.o" "0"
         )
        ADD_TEST(NAME libmmgs_fortran_io_1   COMMAND ${LIBMMGS_EXECFORTRAN_IO}
          "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
          "${CTEST_OUTPUT_DIR}/libmmgs_Fortran_io-torus.o" "1"
         )

      ENDIF()

    ENDIF()

    IF ( ONLY_VERY_SHORT_TESTS )
      SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )

      ADD_TEST(NAME mmgs_very_short   COMMAND ${EXECUT_MMGS}
        "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/cube.mesh"
        "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_0_a-cube.o"
        )

    ELSE ( )
      # Add mmgs tests
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmgs_tests.cmake )
      INCLUDE( ${PROJECT_SOURCE_DIR}/cmake/testing/mmg_tests.cmake )
    ENDIF ( )

  ENDIF ( MMGS_CI )

ENDIF ( BUILD_TESTING )
