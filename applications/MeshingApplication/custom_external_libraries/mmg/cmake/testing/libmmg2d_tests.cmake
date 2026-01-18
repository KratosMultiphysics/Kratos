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
#####        Build executables
#####
###############################################################################

SET ( MMG2D_LIB_TESTS
  libmmg2d_adp_example0_a
  libmmg2d_adp_example0_b
  libmmg2d_adp_example1
  libmmg2d_adp_example2
  libmmg2d_gene_example0
  libmmg2d_ls_example0
  libmmg2d_lsOnly
  libmmg2d_lsOnly_optim
  libmmg2d_lsOnly_hsiz
  libmmg2d_lsAndMetric_optim
  libmmg2d_lsAndMetric_hsiz
  libmmg2d_lsAndMetric
  )

# Additional tests that needs to download ci meshes
IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
  LIST ( APPEND MMG2D_LIB_TESTS test_api2d_0 )
ENDIF ( )

SET ( MMG2D_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_b/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example2/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/squareGeneration_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main_optim.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main_hsiz.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main_optim.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main_hsiz.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main.c
  )

# Additional tests that needs to download ci meshes
#
# Remark: as they are piece of code, it would be probably better to keep it in
# the mmg repository (so we have versionning). Other pieces of code
# (compare-para-tria and ridge-preservation-in-ls-mode) are hosted in the repo,
# I don't think that there is a reason for the difference in the choice of
# hosting...

IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
  LIST ( APPEND MMG2D_LIB_TESTS_MAIN_PATH
    ${MMG2D_CI_TESTS}/API_tests/2d.c
    )
ENDIF( )

IF ( LIBMMG2D_STATIC )
  SET ( lib_name lib${PROJECT_NAME}2d_a )
  SET ( lib_type "STATIC" )
ELSEIF ( LIBMMG2D_SHARED )
  SET ( lib_name lib${PROJECT_NAME}2d_so )
  SET ( lib_type "SHARED" )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)
  ENABLE_LANGUAGE ( Fortran )

  SET ( MMG2D_LIB_TESTS ${MMG2D_LIB_TESTS}
    libmmg2d_fortran_a
    libmmg2d_fortran_b
    libmmg2d_fortran_io
    libmmg2d_fortran_lsOnly
    libmmg2d_fortran_lsOnly_optim
    libmmg2d_fortran_lsOnly_hsiz
    libmmg2d_fortran_lsAndMetric_optim
    libmmg2d_fortran_lsAndMetric_hsiz
    libmmg2d_fortran_lsAndMetric
    test_io2d_fortran
    )
  # Additional tests that needs to download ci meshes
  IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
    LIST ( APPEND MMG2D_LIB_TESTS test_api2d_fortran_0 )
  ENDIF( )

  SET ( MMG2D_LIB_TESTS_MAIN_PATH ${MMG2D_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_b/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main_optim.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/main_hsiz.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main_optim.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main_hsiz.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsAndMetric/main.F90
    ${PROJECT_SOURCE_DIR}/cmake/testing/code/mmg2d_io.F90
    )
  # Additional tests that needs to download ci meshes
  IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
    LIST ( APPEND MMG2D_LIB_TESTS_MAIN_PATH
      ${MMG2D_CI_TESTS}/API_tests/2d.F90
      )
  ENDIF( )


ENDIF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)

LIST(LENGTH MMG2D_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMG2D_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMG2D_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_2d_headers ${lib_name} ${lib_type})

ENDFOREACH ( )

SET ( src_test_met2d
  ${PROJECT_SOURCE_DIR}/src/common/bezier.c
  ${PROJECT_SOURCE_DIR}/src/common/eigenv.c
  ${PROJECT_SOURCE_DIR}/src/common/mettools.c
  ${PROJECT_SOURCE_DIR}/src/common/anisosiz.c
  ${PROJECT_SOURCE_DIR}/src/common/isosiz.c
  ${PROJECT_SOURCE_DIR}/src/common/tools.c
  ${PROJECT_SOURCE_DIR}/src/common/mmgexterns.c
  ${PROJECT_SOURCE_DIR}/cmake/testing/code/test_met2d.c
  )
ADD_LIBRARY_TEST ( test_met2d "${src_test_met2d}" copy_2d_headers ${lib_name} ${lib_type})
TARGET_LINK_LIBRARIES ( test_met2d PRIVATE ${M_LIB} )

###############################################################################
#####
#####  Add tests if needed
#####
###############################################################################

ADD_TEST(NAME libmmg2d_adp_example0_a
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example0_a
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/init.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_a-init.o"
  )
ADD_TEST(NAME libmmg2d_adp_example0_b
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example0_b
  "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_0_b.o.mesh"
  )
ADD_TEST(NAME libmmg2d_adp_example1
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example1
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/dom.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_1-dom.o"
  )
ADD_TEST(NAME libmmg2d_adp_example2
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_adp_example2
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/dom.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_2-dom.o"
  "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_2-dom-end.o"
  )
ADD_TEST(NAME libmmg2d_gene_example0
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_gene_example0
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/squareGeneration_example0/carretest.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_Generation_0-carre.o"
  )
ADD_TEST(NAME libmmg2d_ls0_io_0
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_ls_example0
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_io_0-naca.o" "0"
  )
ADD_TEST(NAME libmmg2d_ls0_io_1
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_ls_example0
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
  "${CTEST_OUTPUT_DIR}/libmmg2d_io_0-naca.o" "1"
  )
ADD_TEST(NAME libmmg2d_lsOnly
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_lsOnly
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
  "${CTEST_OUTPUT_DIR}/libmmg2d_lsOnly_multimat.o"
  )
ADD_TEST(NAME libmmg2d_lsAndMetric
  COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_lsAndMetric
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
  "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
  "${CTEST_OUTPUT_DIR}/libmmg2d_lsAndMetric_multimat.o"
  )

IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
  ADD_TEST(NAME test_api2d_0
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_api2d_0
    "${MMG2D_CI_TESTS}/API_tests/2dom.mesh"
    "${CTEST_OUTPUT_DIR}/test_API2d.o"
    )
ENDIF()

ADD_TEST(NAME test_met2d   COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_met2d  )

IF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)
  ADD_TEST(NAME libmmg2d_fortran_a
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_a
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/init.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d-Adaptation_Fortran_0_a-init.o"
    )
  ADD_TEST(NAME libmmg2d_fortran_b
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_b
    "${CTEST_OUTPUT_DIR}/libmmg2d_Adaptation_Fortran_0_b.o"
    )
  ADD_TEST(NAME libmmg2d_fortran_io_0
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_io
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Fortran_io-naca.o" "0"
    )
  ADD_TEST(NAME libmmg2d_fortran_io_1
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_io
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example0/naca-multiSols.mesh"
    "${CTEST_OUTPUT_DIR}/libmmg2d_Fortran_io-naca.o" "1"
    )
  ADD_TEST(NAME libmmg2d_fortran_lsOnly
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_lsOnly
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmg2d_lsOnly_multimat.o"
    )
  ADD_TEST(NAME libmmg2d_fortran_lsAndMetric
    COMMAND ${EXECUTABLE_OUTPUT_PATH}/libmmg2d_fortran_lsAndMetric
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat.mesh"
    "${PROJECT_SOURCE_DIR}/libexamples/mmg2d/IsosurfDiscretization_lsOnly/multi-mat-sol.sol"
    "${CTEST_OUTPUT_DIR}/libmmg2d_lsAndMetric_multimat.o"
    )

  IF ( MMG2D_CI AND NOT ONLY_VERY_SHORT_TESTS )
    ADD_TEST(NAME test_api2d_fortran_0
      COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_api2d_fortran_0
      "${MMG2D_CI_TESTS}/API_tests/2dom.mesh"
      "${CTEST_OUTPUT_DIR}/test_API2d.o"
      )
    ADD_TEST(NAME test_io2d_fortran_scalar
      COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_io2d_fortran
      "${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh"
      "${CTEST_OUTPUT_DIR}/hybrid-2d-scal.o" 0
      )
    ADD_TEST(NAME test_io2d_fortran_array
      COMMAND ${EXECUTABLE_OUTPUT_PATH}/test_io2d_fortran
      "${MMG2D_CI_TESTS}/Hybrid/hybrid.mesh"
      "${CTEST_OUTPUT_DIR}/hybrid-2d-array.o" 1
      )
  ENDIF()

ENDIF ( )
