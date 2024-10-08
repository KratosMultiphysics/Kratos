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
#####  Build executable for Mmgs Library Examples and add tests if needed
#####
###############################################################################

SET ( MMGS_LIB_TESTS
  libmmgs_example0_a
  libmmgs_example0_b
  libmmgs_example1
  libmmgs_example2
  libmmgs_example3
  libmmgs_lsOnly
  libmmgs_lsOnly_optim
  libmmgs_lsOnly_hsiz
  libmmgs_lsAndMetric_optim
  libmmgs_lsAndMetric_hsiz
  libmmgs_lsAndMetric
  )
SET ( MMGS_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_b/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main_optim.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main_hsiz.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main_optim.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main_hsiz.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main.c
  )

IF ( LIBMMGS_STATIC )
  SET ( lib_name lib${PROJECT_NAME}s_a )
  SET ( lib_type "STATIC" )
ELSEIF ( LIBMMGS_SHARED )
  SET ( lib_name lib${PROJECT_NAME}s_so )
  SET ( lib_type "SHARED" )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)
  ENABLE_LANGUAGE ( Fortran )

  SET ( MMGS_LIB_TESTS ${MMGS_LIB_TESTS}
    libmmgs_fortran_a
    libmmgs_fortran_b
    libmmgs_fortran_io
    libmmgs_fortran_lsOnly
    libmmgs_fortran_lsOnly_optim
    libmmgs_fortran_lsOnly_hsiz
    libmmgs_fortran_lsAndMetric_optim
    libmmgs_fortran_lsAndMetric_hsiz
    libmmgs_fortran_lsAndMetric
    )

  SET ( MMGS_LIB_TESTS_MAIN_PATH ${MMGS_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_a/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_b/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main_optim.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/main_hsiz.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main_optim.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main_hsiz.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsAndMetric/main.F90
    )

ENDIF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)

LIST(LENGTH MMGS_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMGS_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMGS_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_s_headers ${lib_name} ${lib_type})

ENDFOREACH ( )

ADD_TEST(NAME libmmgs_example0_a
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_a
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/cube.mesh"
  "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_0_a-cube.o"
  )
ADD_TEST(NAME libmmgs_example0_b
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example0_b
  "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_0_b.o.mesh"
  )
ADD_TEST(NAME libmmgs_example1
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example1
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/2spheres"
  "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_1-2spheres_1.o"
  "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_1-2spheres_2.o"
  )
ADD_TEST(NAME libmmgs_example2
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example2
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/teapot"
  "${CTEST_OUTPUT_DIR}/libmmgs-IsosurfDiscretization_0-teapot.o"
  )
ADD_TEST(NAME libmmgs_example3_io_0
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example3
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
  "${CTEST_OUTPUT_DIR}/libmmgs_io_3-naca.o" "0"
  )
ADD_TEST(NAME libmmgs_example3_io_1
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_example3
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
  "${CTEST_OUTPUT_DIR}/libmmgs_io_3-naca.o" "1"
  )
ADD_TEST(NAME libmmgs_lsOnly
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_lsOnly
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat.mesh"
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
  "${CTEST_OUTPUT_DIR}/libmmgs_lsOnly_multimat.o"
  )
ADD_TEST(NAME libmmgs_lsAndMetric
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_lsAndMetric
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat.mesh"
  "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
  "${CTEST_OUTPUT_DIR}/libmmgs_lsAndMetric_multimat.o"
  )

IF ( CMAKE_Fortran_COMPILER AND PERL_FOUND ) 
  SET(LIBMMGS_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_a)
  SET(LIBMMGS_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_b)
  SET(LIBMMGS_EXECFORTRAN_IO ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_io)
  SET(LIBMMGS_EXECFORTRAN_LSONLY ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_lsOnly )
  SET(LIBMMGS_EXECFORTRAN_LSANDMETRIC ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_lsAndMetric )

  ADD_TEST(NAME libmmgs_fortran_a
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_a
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_a/cube.mesh"
    "${CTEST_OUTPUT_DIR}/libmmgs-Adaptation_Fortran_0_a-cube.o"
    )
  ADD_TEST(NAME libmmgs_fortran_b
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_b
    "${CTEST_OUTPUT_DIR}/libmmgs_Adaptation_Fortran_0_b.o"
    )
  ADD_TEST(NAME libmmgs_fortran_io_0
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_io
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
    "${CTEST_OUTPUT_DIR}/libmmgs_Fortran_io-torus.o" "0"
    )
  ADD_TEST(NAME libmmgs_fortran_io_1
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_io
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/torus.mesh"
    "${CTEST_OUTPUT_DIR}/libmmgs_Fortran_io-torus.o" "1"
    )
  ADD_TEST(NAME libmmgs_fortran_lsOnly
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_lsOnly
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmgs_lsOnly_multimat.o"
    )
  ADD_TEST(NAME libmmgs_fortran_lsAndMetric
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmgs_fortran_lsAndMetric
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmgs_lsAndMetric_multimat.o"
    )

ENDIF ( )
