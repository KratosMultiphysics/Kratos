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
#####         Mmgs Library Examples
#####
###############################################################################

SET ( MMGS_LIB_TESTS
  libmmgs_example0_a
  libmmgs_example0_b
  libmmgs_example1
  libmmgs_example2
  libmmgs_example3
  )
SET ( MMGS_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_a/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0/example0_b/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example1/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/IsosurfDiscretization_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/main.c
  )

IF ( LIBMMGS_STATIC )
  SET ( lib_name lib${PROJECT_NAME}s_a )
ELSEIF ( LIBMMGS_SHARED )
  SET ( lib_name lib${PROJECT_NAME}s_so )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF ( CMAKE_Fortran_COMPILER )
  ENABLE_LANGUAGE ( Fortran )

  SET ( MMGS_LIB_TESTS ${MMGS_LIB_TESTS}
    libmmgs_fortran_a
    libmmgs_fortran_b
    libmmgs_fortran_io
    )

  SET ( MMGS_LIB_TESTS_MAIN_PATH ${MMGS_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_a/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/adaptation_example0_fortran/example0_b/main.F90
    ${PROJECT_SOURCE_DIR}/libexamples/mmgs/io_multisols_example3/main.F90
    )

ENDIF ( CMAKE_Fortran_COMPILER )

LIST(LENGTH MMGS_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMGS_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMGS_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_s_headers ${lib_name} )

ENDFOREACH ( )
