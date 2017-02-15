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

SET(MMG3D_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/mmg3d)

############################################################################
#####
#####         Fortran header: libmmg3df.h
#####
############################################################################

IF ( NOT WIN32 )
  ADD_CUSTOM_COMMAND(OUTPUT ${MMG3D_SOURCE_DIR}/libmmg3df.h
    COMMAND genheader ${MMG3D_SOURCE_DIR}/libmmg3df.h
    ${MMG3D_SOURCE_DIR}/libmmg3d.h ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS genheader ${MMG3D_SOURCE_DIR}/libmmg3d.h
    ${COMMON_SOURCE_DIR}/libmmgtypesf.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for mmg3d"
    )
ENDIF()

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

# Header files
INCLUDE_DIRECTORIES(${MMG3D_SOURCE_DIR})

# Source files
FILE(
  GLOB
  source_files
  ${MMG3D_SOURCE_DIR}/*.c   ${MMG3D_SOURCE_DIR}/*.h
  ${COMMON_SOURCE_DIR}/*.c ${COMMON_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM source_files
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.c
  ${CMAKE_SOURCE_DIR}/src/libmmg.h
  ${CMAKE_SOURCE_DIR}/src/libmmgf.h
  ${REMOVE_FILE})
FILE(
  GLOB
  main_file
  ${MMG3D_SOURCE_DIR}/mmg3d.c
  )
FILE(
  GLOB
  lib_file
  #${MMG3D_SOURCE_DIR}/library_tools.c
  ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.c
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
    "problem (see https://github.com/SUscTools/Elas to download it). "
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
  ADD_LIBRARY(${PROJECT_NAME}3d_a  STATIC
    ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
    ${source_files} ${lib_file} )
  SET_TARGET_PROPERTIES(${PROJECT_NAME}3d_a PROPERTIES OUTPUT_NAME
    ${PROJECT_NAME}3d)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d_a ${LIBRARIES})
  INSTALL(TARGETS ${PROJECT_NAME}3d_a
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

# Compile shared library
IF ( LIBMMG3D_SHARED )
  ADD_LIBRARY(${PROJECT_NAME}3d_so SHARED
    ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
    ${source_files} ${lib_file})
  SET_TARGET_PROPERTIES(${PROJECT_NAME}3d_so PROPERTIES
    VERSION ${CMAKE_RELEASE_VERSION} SOVERSION 5)
  SET_TARGET_PROPERTIES(${PROJECT_NAME}3d_so PROPERTIES
    OUTPUT_NAME ${PROJECT_NAME}3d)
  TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d_so ${LIBRARIES})
  INSTALL(TARGETS ${PROJECT_NAME}3d_so
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib)
ENDIF()

IF ( LIBMMG3D_STATIC OR LIBMMG3D_SHARED )
  # mmg3d header files needed for library
  SET( mmg3d_headers
    ${MMG3D_SOURCE_DIR}/libmmg3d.h
    ${MMG3D_SOURCE_DIR}/libmmg3df.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${COMMON_SOURCE_DIR}/libmmgtypesf.h
    )
  SET(MMG3D_INCLUDE ${CMAKE_SOURCE_DIR}/include/mmg/mmg3d )
  SET( mmg3d_includes
    ${MMG3D_INCLUDE}/libmmg3d.h
    ${MMG3D_INCLUDE}/libmmg3df.h
    ${MMG3D_INCLUDE}/libmmgtypes.h
    ${MMG3D_INCLUDE}/libmmgtypesf.h
    )
  # Install header files in /usr/local or equivalent
  INSTALL(FILES ${mmg3d_headers} DESTINATION include/mmg/mmg3d)

  ADD_CUSTOM_COMMAND(OUTPUT ${MMG3D_INCLUDE}/libmmgtypesf.h
    COMMAND ${CMAKE_COMMAND} -E copy ${COMMON_SOURCE_DIR}/libmmgtypesf.h
    ${MMG3D_INCLUDE}/libmmgtypesf.h
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${COMMON_SOURCE_DIR}/libmmgtypesf.h)
  ADD_CUSTOM_COMMAND(OUTPUT ${MMG3D_INCLUDE}/libmmg3df.h
    COMMAND ${CMAKE_COMMAND} -E copy ${MMG3D_SOURCE_DIR}/libmmg3df.h ${MMG3D_INCLUDE}/libmmg3df.h
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS ${MMG3D_SOURCE_DIR}/libmmg3df.h)

  # Install header files in project directory
  FILE(INSTALL  ${mmg3d_headers} DESTINATION ${MMG3D_INCLUDE}
    PATTERN "libmmg*f.h"  EXCLUDE)

  ADD_CUSTOM_TARGET(copy_3d_headers ALL
    DEPENDS  ${MMG3D_INCLUDE}/libmmg3df.h  ${MMG3D_INCLUDE}/libmmg3d.h
    ${MMG3D_INCLUDE}/libmmgtypesf.h ${MMG3D_INCLUDE}/libmmgtypes.h )

ENDIF()

############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG3D )
  INCLUDE(cmake/testing/libmmg3d_tests.cmake)
  INCLUDE(cmake/testing/libmmg3d_oldAPI_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile MMG3D executable
#####
###############################################################################
ADD_EXECUTABLE(${PROJECT_NAME}3d
  ${MMG3D_SOURCE_DIR}/lib${PROJECT_NAME}3df.h
  ${source_files} ${main_file})

IF ( WIN32 AND NOT MINGW AND USE_SCOTCH )
  my_add_link_flags(${PROJECT_NAME}3d "/SAFESEH:NO")
ENDIF ( )

TARGET_LINK_LIBRARIES(${PROJECT_NAME}3d ${LIBRARIES})
INSTALL(TARGETS ${PROJECT_NAME}3d RUNTIME DESTINATION bin)

# in debug mode we name the executable mmg3d_debug
SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES DEBUG_POSTFIX _debug)
# in Release mode we name the executable mmg3d_O3
SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES RELEASE_POSTFIX _O3)
# in RelWithDebInfo mode we name the executable mmg3d_O3d
SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES RELWITHDEBINFO_POSTFIX _O3d)
# in MinSizeRel mode we name the executable mmg3d_Os
SET_TARGET_PROPERTIES(${PROJECT_NAME}3d PROPERTIES MINSIZEREL_POSTFIX _Os)


###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

IF ( BUILD_TESTING )
  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET(MMG3D_CI_TESTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg3d )
  SET(MMG_CI_TESTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg )

  ##-------------------------------------------------------------------##
  ##--------------------------- Add tests and configure it ------------##
  ##-------------------------------------------------------------------##
  # Add runtime that we want to test for mmg3d
  IF ( MMG3D_CI )

    IF ( LONG_TESTS )
      # Run some tests twice with the output of the previous test as input
      OPTION ( RUN_AGAIN "Enable/Disable second run of some tests" ON )
      MARK_AS_ADVANCED(RUN_AGAIN)
    ELSE ( )
      SET ( RUN_AGAIN OFF )
    ENDIF ( )

    IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
      SET(EXECUT_MMG3D ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_debug)
      SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
    ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
      SET(EXECUT_MMG3D ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_O3)
      SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
    ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
      SET(EXECUT_MMG3D ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_O3d)
      SET(BUILDNAME ${BUILDNAME}_O3d CACHE STRING "build name variable")
    ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "MinSizeRel")
      SET(EXECUT_MMG3D ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d_Os)
      SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
    ELSE()
      SET(EXECUT_MMG3D ${EXECUTABLE_OUTPUT_PATH}/${PROJECT_NAME}3d)
      SET(BUILDNAME ${BUILDNAME} CACHE STRING "build name variable")
    ENDIF()

    SET ( LISTEXEC_MMG3D ${EXECUT_MMG3D} )
    IF ( TEST_LIBMMG3D )
      SET(LIBMMG3D_EXEC0_a ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_a)
      SET(LIBMMG3D_EXEC0_b ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example0_b)
      SET(LIBMMG3D_EXEC1   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example1)
      SET(LIBMMG3D_EXEC2   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example2)
      SET(LIBMMG3D_EXEC4   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example4)
      SET(LIBMMG3D_EXEC5   ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_example5)

      ADD_TEST(NAME libmmg3d_example0_a COMMAND ${LIBMMG3D_EXEC0_a})
      ADD_TEST(NAME libmmg3d_example0_b COMMAND ${LIBMMG3D_EXEC0_b})
      ADD_TEST(NAME libmmg3d_example1   COMMAND ${LIBMMG3D_EXEC1})
      ADD_TEST(NAME libmmg3d_example2   COMMAND ${LIBMMG3D_EXEC2})
      IF ( USE_ELAS )
        ADD_TEST(NAME libmmg3d_example4   COMMAND ${LIBMMG3D_EXEC4})
      ENDIF ()
      ADD_TEST(NAME libmmg3d_example5   COMMAND ${LIBMMG3D_EXEC5})

      SET( LISTEXEC_MMG3D ${LISTEXEC_MMG3D} )

      IF ( CMAKE_Fortran_COMPILER)
        SET(LIBMMG3D_EXECFORTRAN_a ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_a)
        SET(LIBMMG3D_EXECFORTRAN_b ${EXECUTABLE_OUTPUT_PATH}/libmmg3d_fortran_b)
        ADD_TEST(NAME libmmg3d_fortran_a  COMMAND ${LIBMMG3D_EXECFORTRAN_a})
        ADD_TEST(NAME libmmg3d_fortran_b  COMMAND ${LIBMMG3D_EXECFORTRAN_b})
      ENDIF()

    ENDIF ( TEST_LIBMMG3D )

    # Add more tests
    INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg3d_tests.cmake )

    IF ( RUN_AGAIN )
      INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg3d_rerun_tests.cmake )
    ENDIF()

  ENDIF ( MMG3D_CI )

ENDIF ( BUILD_TESTING )
