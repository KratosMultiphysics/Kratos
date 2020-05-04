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
#####        Mmg2d library examples
#####
###############################################################################

SET ( MMG2D_LIB_TESTS
  libmmg2d_example0_a
  libmmg2d_example0_b
  libmmg2d_example1
  libmmg2d_example2
  libmmg2d_example3
  test_api2d_0
  )
SET ( MMG2D_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_a/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0/example0_b/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example1/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/squareGeneration_example2/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example3/main.c
  ${MMG2D_CI_TESTS}/API_tests/2d.c
  )

IF ( LIBMMG2D_STATIC )
  SET ( lib_name lib${PROJECT_NAME}2d_a )
ELSEIF ( LIBMMG2D_SHARED )
  SET ( lib_name lib${PROJECT_NAME}2d_so )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF ( CMAKE_Fortran_COMPILER )
  ENABLE_LANGUAGE ( Fortran )

  SET ( MMG2D_LIB_TESTS ${MMG2D_LIB_TESTS}
    libmmg2d_fortran_a
    libmmg2d_fortran_b
    libmmg2d_fortran_io
    test_api2d_fortran_0
    )

  SET ( MMG2D_LIB_TESTS_MAIN_PATH ${MMG2D_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_a/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/adaptation_example0_fortran/example0_b/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmg2d/io_multisols_example3/main.F90
    ${MMG2D_CI_TESTS}/API_tests/2d.F90
    )

ENDIF ( CMAKE_Fortran_COMPILER )

LIST(LENGTH MMG2D_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMG2D_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMG2D_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_2d_headers ${lib_name} )

ENDFOREACH ( )
