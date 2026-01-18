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
#####         Mmg Library Examples
#####
###############################################################################

SET ( MMG_LIB_TESTS
  libmmg_example0_a
  libmmg_cpp_a
  )

SET ( MMG_LIB_TESTS_MAIN_PATH
  ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0/main.c
  ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_cpp/main.cpp
  )

IF ( LIBMMG_STATIC )
  SET ( lib_name lib${PROJECT_NAME}_a )
  SET ( lib_type "STATIC" )
ELSEIF ( LIBMMG_SHARED )
  SET ( lib_name lib${PROJECT_NAME}_so )
  SET ( lib_type "SHARED" )
ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ( )

#####         Fortran Tests
IF (CMAKE_Fortran_COMPILER AND PERL_FOUND)
  ENABLE_LANGUAGE (Fortran)

  SET ( MMG_LIB_TESTS ${MMG_LIB_TESTS}
    libmmg_fortran_a
    )

  SET ( MMG_LIB_TESTS_MAIN_PATH ${MMG_LIB_TESTS_MAIN_PATH}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/main.F90
    )

  ENDIF (CMAKE_Fortran_COMPILER AND PERL_FOUND)

LIST(LENGTH MMG_LIB_TESTS nbTests_tmp)
MATH(EXPR nbTests "${nbTests_tmp} - 1")

FOREACH ( test_idx RANGE ${nbTests} )
  LIST ( GET MMG_LIB_TESTS           ${test_idx} test_name )
  LIST ( GET MMG_LIB_TESTS_MAIN_PATH ${test_idx} main_path )

  ADD_LIBRARY_TEST ( ${test_name} ${main_path} copy_mmg_headers ${lib_name} ${lib_type} )

ENDFOREACH ( )

IF( MMG_CI )
  # Add libmmg tests
  SET(LIBMMG_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_example0_a )
  SET(LIBMMG_CPP_a   ${EXECUTABLE_OUTPUT_PATH}/libmmg_cpp_a )

  ADD_TEST(NAME libmmg_example0_a   COMMAND ${LIBMMG_EXEC0_a}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0/init
    ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0/cube
    "${CTEST_OUTPUT_DIR}/libmmg_Adaptation_0.o")
  ADD_TEST(NAME libmmg_cpp_a        COMMAND ${LIBMMG_CPP_a}
    ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_cpp/init
    ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_cpp/cube
    "${CTEST_OUTPUT_DIR}/libmmg_Adaptation_0_cpp.o")

  IF ( CMAKE_Fortran_COMPILER AND PERL_FOUND)
    SET(LIBMMG_FORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_fortran_a)
    ADD_TEST(NAME libmmg_fortran   COMMAND ${LIBMMG_FORTRAN_a}
      ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/init
      ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/cube
      "${CTEST_OUTPUT_DIR}/libmmg_Adaptation_0_Fortran.o"
      )
  ENDIF()

ENDIF( MMG_CI )
