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
  ${MMG2D_SOURCE_DIR}/*.c
  ${MMG3D_SOURCE_DIR}/*.c
  ${MMGS_SOURCE_DIR}/*.c
  ${COMMON_SOURCE_DIR}/*.c
  )
LIST(REMOVE_ITEM mmg_library_files
  ${MMG2D_SOURCE_DIR}/mmg2d.c
  ${MMGS_SOURCE_DIR}/mmgs.c
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  ${REMOVE_FILE} )


IF ( LIBMMG_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_a  STATIC
    "${mmg_library_files}"
    ${PROJECT_NAME} )
ENDIF()

# Compile shared library
IF ( LIBMMG_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_so  SHARED
    "${mmg_library_files}"
    ${PROJECT_NAME} )
ENDIF()

IF ( LIBMMG_STATIC OR LIBMMG_SHARED )
  # mmg header files needed for library
  SET( mmg2d_headers
    ${MMG2D_SOURCE_DIR}/libmmg2d.h
    ${MMG2D_BINARY_DIR}/libmmg2df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmg3d_headers
    ${MMG3D_SOURCE_DIR}/libmmg3d.h
    ${MMG3D_BINARY_DIR}/libmmg3df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmgs_headers
    ${MMGS_SOURCE_DIR}/libmmgs.h
    ${MMGS_BINARY_DIR}/libmmgsf.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_BINARY_DIR}/libmmgtypesf.h
    )
  SET( mmg_headers
    ${PROJECT_SOURCE_DIR}/src/mmg/libmmg.h
    ${PROJECT_SOURCE_DIR}/src/mmg/libmmgf.h
    )
  SET(MMG2D_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmg2d )
  SET(MMGS_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmgs )
  SET(MMG3D_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg/mmg3d )
  SET(MMG_INCLUDE ${PROJECT_BINARY_DIR}/include/mmg )

  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg2d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg2d)
  INSTALL(FILES ${mmgs_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmgs)
  INSTALL(FILES ${mmg3d_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg/mmg3d)
  INSTALL(FILES ${mmg_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mmg)

  FILE(INSTALL ${PROJECT_SOURCE_DIR}/src/mmg/libmmgf.h DESTINATION  ${PROJECT_BINARY_DIR}/include/mmg/)

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


  ADD_CUSTOM_TARGET(copy_mmg_headers ALL
    DEPENDS
    copy_2d_headers copy_s_headers copy_3d_headers
    ${PROJECT_BINARY_DIR}/include/mmg/libmmgf.h
    ${PROJECT_BINARY_DIR}/include/mmg/libmmg.h
    ${PROJECT_BINARY_DIR}/include/mmg/mmg3d/libmmgtypes.h )

ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG )
  INCLUDE(cmake/testing/libmmg_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )

  SET ( CTEST_OUTPUT_DIR ${PROJECT_BINARY_DIR}/TEST_OUTPUTS )
  FILE ( MAKE_DIRECTORY  ${CTEST_OUTPUT_DIR} )

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg
  IF( MMG_CI )
    # Add libmmg tests
    IF ( TEST_LIBMMG )
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

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG_FORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg_fortran_a)
        ADD_TEST(NAME libmmg_fortran   COMMAND ${LIBMMG_FORTRAN_a}
          ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/init
          ${PROJECT_SOURCE_DIR}/libexamples/mmg/adaptation_example0_fortran/cube
          "${CTEST_OUTPUT_DIR}/libmmg_Adaptation_0_Fortran.o"
          )
      ENDIF()
    ENDIF ()

  ENDIF( MMG_CI )

ENDIF ( BUILD_TESTING )
