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

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.0)
INCLUDE(CMakeDependentOption)

PROJECT (mmg)

###############################################################################
#####
#####         Release version and date
#####
###############################################################################
SET (CMAKE_RELEASE_VERSION_MAJOR "5" )
SET (CMAKE_RELEASE_VERSION_MINOR "2" )
SET (CMAKE_RELEASE_VERSION_PATCH "2" )
SET (CMAKE_RELEASE_DATE "Feb 8, 2017" )

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

###############################################################################
#####
#####         Set preprocessor flag to detect windows architecture and MinGW
#####
###############################################################################
IF (NOT WIN32)
  # Set preprocessor flags to say that we are posix and gnu compatible
  SET (DEF_POSIX "#define POSIX")
  SET (DEF_GNU "#define GNU")
ELSEIF(MINGW)
  SET (DEF_GNU "#define GNU")
ENDIF ()

# Create mmgcommon.h file with the good release and architecture infos.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgcommon.h.in
  ${CMAKE_CURRENT_SOURCE_DIR}/src/common/mmgcommon.h @ONLY)

###############################################################################
#####
#####         Main CMake compilation variables
#####
###############################################################################

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
SET(COMMON_SOURCE_DIR ${CMAKE_SOURCE_DIR}/src/common)

# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
  # using clang
  SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
  # using icc
  SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()

# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS " -Wno-char-subscripts ${CMAKE_C_FLAGS}")
  IF(APPLE)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64 ${CMAKE_EXE_LINKER_FLAGS}")

    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
    #  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF()
ENDIF()

CMAKE_DEPENDENT_OPTION ( PROFILING "Enable/Disable PROFILING" OFF
  "CMAKE_COMPILER_IS_GNUCC" OFF )

IF(PROFILING)
  ADD_DEFINITIONS(-pg)
  SET(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
ENDIF(PROFILING)

###############################################################################
#####
#####         Choose executable target to compile
#####
###############################################################################
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  seT(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
  SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE})
ENDIF()

SET ( BUILD "MMG" CACHE STRING
  "\"MMG3D\", \"MMGS\",\"MMG2D\", or \"MMG\" to compile both MMGS, MMG3D and MMG2D" )

IF ( BUILD MATCHES "MMGS" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       ON  )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG3D" )
  SET(BUILD_MMG3D      ON  )
  SET(BUILD_MMGS       OFF )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG2D      OFF  )
ELSEIF ( BUILD MATCHES "MMG2D" )
  SET(BUILD_MMG        OFF )
  SET(BUILD_MMG3D      OFF )
  SET(BUILD_MMGS       OFF  )
  SET(BUILD_MMG2D      ON  )
ELSEIF ( BUILD MATCHES "MMG" )
  SET(BUILD_MMG        ON  )
  SET(BUILD_MMGS       ON )
  SET(BUILD_MMG3D      ON )
  SET(BUILD_MMG2D      ON  )
ELSE ( )
  MESSAGE ( ERROR " Nothing to build.\n"
    "Please choose the target to build (set the BUILD CMake's"
    " variable to \"MMG\", \"MMG2D\", \"MMGS\" OR"
    " \"MMG3D\")" )

ENDIF()

IF ( BUILD_MMG3D OR BUILD_MMGS )
 SET(BUILD_MMGS3D ON)
ELSE ()
  SET(BUILD_MMGS3D OFF)
ENDIF()
############################################################################
#####
#####         Fortran header: libmmgcommonf.h
#####
############################################################################

IF ( NOT WIN32 )
  ADD_EXECUTABLE(genheader ${CMAKE_SOURCE_DIR}/scripts/genheader.c)
  ADD_CUSTOM_COMMAND(OUTPUT ${COMMON_SOURCE_DIR}/libmmgtypesf.h
    COMMAND genheader ${COMMON_SOURCE_DIR}/libmmgtypesf.h
    ${COMMON_SOURCE_DIR}/libmmgtypes.h ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    DEPENDS genheader ${COMMON_SOURCE_DIR}/libmmgtypes.h
    ${CMAKE_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for libmmgtypes.h file"
    )
ENDIF()

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${COMMON_SOURCE_DIR})

# Libraries
IF ( NOT WIN32 OR MINGW )
  FIND_LIBRARY(M_LIB m)
  MARK_AS_ADVANCED(FORCE M_LIB)
  SET( LIBRARIES ${M_LIB} ${LIBRARIES})
ENDIF()

############################################################################
#####
#####         Scotch
#####
############################################################################
# add SCOTCH library?
IF( USE_SCOTCH )
  # Set flags for building test program
  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIR})

  SET(CMAKE_REQUIRED_INCLUDES ${SCOTCH_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCH_LIBRARY})
  SET(CMAKE_REQUIRED_LIBRARIES ${SCOTCHERR_LIBRARY})

  SET(CMAKE_C_FLAGS "-DUSE_SCOTCH ${CMAKE_C_FLAGS}")

  IF ( NOT WIN32 )
    SET(SCOTCH_LINK_FLAGS "-lpthread")
  ELSE ()
    # On win32, it is easier to compile scotch without the pthreads so, don't use the -lpthread flag

    # Turn off "Image Has Safe Exception Handlers" linker option (SAFESEH error of visual studio):
    # see the my_add_link_flags macro.
    INCLUDE(cmake/config/AddCompileLinkFlags.cmake)
  ENDIF ()

  MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARY}")
  SET( LIBRARIES ${SCOTCH_LINK_FLAGS} ${LIBRARIES})
  SET( LIBRARIES ${SCOTCHERR_LIBRARY} ${LIBRARIES})
  SET( LIBRARIES ${SCOTCH_LIBRARY} ${LIBRARIES})
ENDIF()

IF ( SCOTCH_NOTFOUND )
  MESSAGE ( WARNING "Scotch library not found:"
    "Using scotch reduce the execution time of mmg3d "
    "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
    "If you have already installed Scotch and want to use it, "
    "please set the CMake variable or environment variable SCOTCH_DIR "
    "to your scotch directory.")
ENDIF ( )

############################################################################
#####
##### RPATH for MacOSX
#####
############################################################################

SET(CMAKE_MACOSX_RPATH 1)

###############################################################################
#####
#####         Add dependent options
#####
###############################################################################

#------------------------------- mmg2d
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_STATIC "Compile the mmg2d static library" ON
  "BUILD_MMG2D" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG2D_SHARED "Compile the mmg2d dynamic library" OFF
  "BUILD_MMG2D" OFF )
IF (LIBMMG2D_STATIC OR LIBMMG2D_SHARED)
  SET(LIBMMG2D_INTERNAL ON )
ELSE()
  SET(LIBMMG2D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG2D "Compile tests for the mmg2d library" OFF
  "BUILD_MMG2D;LIBMMG2D_INTERNAL" OFF
  )

CMAKE_DEPENDENT_OPTION(
  MMG2D_CI "Enable/Disable continuous integration for mmg2d" ON
  "BUILD_TESTING;BUILD_MMG2D" OFF)

#------------------------------- mmgs
CMAKE_DEPENDENT_OPTION ( LIBMMGS_STATIC "Compile the mmgs static library" ON
  "BUILD_MMGS" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMGS_SHARED "Compile the mmgs dynamic library" OFF
  "BUILD_MMGS" OFF)
IF (LIBMMGS_STATIC OR LIBMMGS_SHARED)
  SET(LIBMMGS_INTERNAL ON )
ELSE()
  SET(LIBMMGS_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMGS "Compile tests for the mmgs library" OFF
  "LIBMMGS_INTERNAL;BUILD_MMGS" OFF
  )

CMAKE_DEPENDENT_OPTION(
  MMGS_CI "Enable/Disable continuous integration for mmgs" ON
  "BUILD_TESTING;BUILD_MMGS" OFF )

#------------------------------- mmg3d
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_STATIC "Compile the mmg3d static library" ON
  "BUILD_MMG3D" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG3D_SHARED "Compile the mmg3d dynamic library" OFF
  "BUILD_MMG3D" OFF )
IF (LIBMMG3D_STATIC OR LIBMMG3D_SHARED)
  SET(LIBMMG3D_INTERNAL ON )
ELSE()
  SET(LIBMMG3D_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG3D "Compile tests for the mmg3d library" OFF
  "BUILD_MMG3D;LIBMMG3D_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG3D_CI "Enable/Disable continuous integration for mmg3d" ON
  "BUILD_TESTING;BUILD_MMG3D" OFF)

CMAKE_DEPENDENT_OPTION ( PATTERN
  "If ON, insertion by patterns, otherwise insertion by delaunay kernel"
  OFF
  "BUILD_MMG3D" OFF)
#-- Remove the next line to have this option visible in basic cmake mode
MARK_AS_ADVANCED(PATTERN)

# add Elastic library?
SET(ELAS_DIR "" CACHE PATH "Installation directory for Elas")
INCLUDE(cmake/modules/FindElas.cmake)
CMAKE_DEPENDENT_OPTION (
  USE_ELAS "Use the Elas library for lagrangian motion option" ON
  "BUILD_MMG3D;ELAS_FOUND" OFF)

# add Scotch library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")
INCLUDE(cmake/modules/FindScotch.cmake)

CMAKE_DEPENDENT_OPTION ( USE_SCOTCH "Use SCOTCH TOOL for renumbering" ON
  "SCOTCH_FOUND;BUILD_MMGS3D" OFF)

#---------------------------- mmg library gathering mmg2d, mmgs and mmg3d libs
CMAKE_DEPENDENT_OPTION ( LIBMMG_STATIC
  "Compile the mmg static library (mmg2d + mmgs + mmg3d)" ON
  "BUILD_MMG" OFF)
CMAKE_DEPENDENT_OPTION ( LIBMMG_SHARED
  "Compile the mmg dynamic library (mmg2d + mmgs + mmg3d)" OFF
  "BUILD_MMG" OFF )
IF (LIBMMG_STATIC OR LIBMMG_SHARED)
  SET(LIBMMG_INTERNAL ON )
ELSE()
  SET(LIBMMG_INTERNAL OFF )
ENDIF()
CMAKE_DEPENDENT_OPTION (
  TEST_LIBMMG "Compile tests for the mmg library" OFF
  "BUILD_MMG;LIBMMG_INTERNAL" OFF
  )
CMAKE_DEPENDENT_OPTION(
  MMG_CI "Enable/Disable continuous integration for mmg" ON
  "BUILD_TESTING;BUILD_MMG" OFF)


###############################################################################
#####
#####         Build Targets and libmmg if needed
#####
###############################################################################
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

IF ( BUILD_MMG2D )
  INCLUDE(cmake/modules/mmg2d.cmake)
ENDIF ( )

IF ( BUILD_MMGS )
  INCLUDE(cmake/modules/mmgs.cmake)
ENDIF ( )

IF ( BUILD_MMG3D )
  INCLUDE(cmake/modules/mmg3d.cmake)
ENDIF ( )

IF ( BUILD_MMG )
  INCLUDE(cmake/modules/mmg.cmake)
ENDIF ( )

###############################################################################
#####
#####         Continuous integration
#####
###############################################################################

OPTION ( BUILD_TESTING "Enable/Disable continuous integration" OFF )

CMAKE_DEPENDENT_OPTION (
  LONG_TESTS "Enable/Disable time consuming tests for the project" OFF
  "BUILD_TESTING" OFF )
MARK_AS_ADVANCED(LONG_TESTS)

CMAKE_DEPENDENT_OPTION ( COVERAGE "Watch code coverage (Linux only)" OFF
"BUILD_TESTING;CMAKE_COMPILER_IS_GNUCC" OFF )

IF( BUILD_TESTING )

  ##-------------------------------------------------------------------##
  ##------- Get the continuous integration tests ----------------------##
  ##-------------------------------------------------------------------##
  INCLUDE(cmake/modules/LoadCiTests.cmake)

  ##-------------------------------------------------------------------##
  ##------- Set the continuous integration options --------------------##
  ##-------------------------------------------------------------------##
  SET(MMG_CI_TESTS ${CMAKE_SOURCE_DIR}/ci_tests/mmg )

  # Remove the --show-reachable=yes option from the default valgrind command
  SET(VALGRIND_COMMAND_OPTIONS "-q --tool=memcheck ")
  SET(VALGRIND_COMMAND_OPTIONS "--leak-check=yes ${VALGRIND_COMMAND_OPTIONS}")
  SET(VALGRIND_COMMAND_OPTIONS "--workaround-gcc296-bugs=yes ${VALGRIND_COMMAND_OPTIONS}")
  SET(VALGRIND_COMMAND_OPTIONS "--num-callers=100 ${VALGRIND_COMMAND_OPTIONS}")

  # Ignore errors coming from the libmesh library
  IF(APPLE)
    # and errors from the libSystem.B library
    SET(MEMORYCHECK_SUPPRESSIONS_FILE
      "${CMAKE_SOURCE_DIR}/cmake/Supp/libSystem.supp")
  ENDIF(APPLE)
  # Copy config file and custom file in the good directories
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/CTestConfig.cmake
    ${CMAKE_SOURCE_DIR}/CTestConfig.cmake)
  CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/cmake/config/CTestCustom.cmake
    ./CTestCustom.cmake)

  # Max time for test case
  SET(DART_TESTING_TIMEOUT 2500 FORCE)

  ENABLE_TESTING()
  INCLUDE(CTest)

  IF ( COVERAGE )
    # Add coverage flags
    SET(CMAKE_C_FLAGS_DEBUG " -fprofile-arcs ${CMAKE_C_FLAGS_DEBUG}")
    SET(CMAKE_C_FLAGS_DEBUG " -ftest-coverage ${CMAKE_C_FLAGS_DEBUG}")
    SET(CMAKE_C_FLAGS_DEBUG " --coverage ${CMAKE_C_FLAGS_DEBUG}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG " -fprofile-arcs ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG " -ftest-coverage ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG " --coverage ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
  ENDIF()

  IF ( CMAKE_COMPILER_IS_GNUCC )
    SET(CMAKE_C_FLAGS_DEBUG "-g -Wall ${CMAKE_C_FLAGS_DEBUG}")
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
  ELSEIF( CMAKE_COMPILER_IS_CLANG )
    SET(CMAKE_C_FLAGS_DEBUG "-g -Weverything -Wno-sign-conversion")
    SET(${CMAKE_C_FLAGS_DEBUG} "-Wno-char-subscripts -Wno-padded "
      ${CMAKE_C_FLAGS_DEBUG})
    SET(CMAKE_EXE_LINKER_FLAGS_DEBUG "-g ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
  ENDIF ( )

  IF ( MMG3D_CI OR MMGS_CI )
    SET(LISTEXEC_MMG ${LISTEXEC_MMGS} ${LISTEXEC_MMG3D})
    SET(EXECUT_MMG ${EXECUT_MMGS} ${EXECUT_MMG3D})

    # Change authorisations for some files...
    IF ( EXISTS ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/unreadable.meshb)
      FILE(REMOVE ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)
    ENDIF()
    FILE(COPY ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/dout.mesh
      DESTINATION ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd7/unwrittable.meshb)

    IF ( EXISTS ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)
      FILE(REMOVE ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)
    ENDIF()
    FILE(COPY ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/dsol.sol
      DESTINATION ${MMG_CI_TESTS}/LeakCheck_AbnormalEnd8/unwrittable.sol)

    # Add common tests for mmgs/3d appli
    INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/mmg_tests.cmake )

  ENDIF(MMG3D_CI OR MMGS_CI)

ENDIF ( BUILD_TESTING )

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)

  # MMG2D Documentation
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/Doxyfile @ONLY)

  ADD_CUSTOM_TARGET(mmg2d_doc
    COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/Doxyfile
    # COMMAND ${CMAKE_COMMAND} -E chdir
    # ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/latex/refman.pdf
    # ${CMAKE_CURRENT_SOURCE_DIR} WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d COMMENT "Generating MMG2D API
    documentation with Doxygen.  Open up the
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg2d/html/index.html file to see
    it." VERBATIM )

  # MMGS Documentation
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile @ONLY)
  ADD_CUSTOM_TARGET(mmgs_doc
    COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/Doxyfile
    # COMMAND ${CMAKE_COMMAND} -E chdir
    # ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/latex/refman.pdf
    # ${CMAKE_CURRENT_SOURCE_DIR} WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs COMMENT "Generating MMGS API
    documentation with Doxygen.  Open up the
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmgs/html/index.html file to see
    it." VERBATIM )

  # MMG3D Documentation
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile @ONLY)
  ADD_CUSTOM_TARGET(mmg3d_doc
    COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/Doxyfile
    # COMMAND ${CMAKE_COMMAND} -E chdir
    # ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/latex/refman.pdf
    # ${CMAKE_CURRENT_SOURCE_DIR} WORKING_DIRECTORY
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d
    COMMENT "Generating MMG3D API documentation with Doxygen. Open up the
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg3d/html/index.html file to see
    it" VERBATIM )

  ADD_CUSTOM_TARGET(doc DEPENDS mmgs_doc mmg3d_doc mmg2d_doc COMMENT "Generating
    MMG API documentation with Doxygen.  Open up the
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/mmg*/html/index.html file to see it"
    VERBATIM )
ENDIF ( DOXYGEN_FOUND )

###############################################################################
#####
#####         Create Package
#####
###############################################################################
INCLUDE(InstallRequiredSystemLibraries)
SET( CPACK_SOURCE_GENERATOR TGZ )
SET( CPACK_GENERATOR TGZ )

SET(CPACK_PACKAGE_NAME ${PROJECT_NAME})
SET(CPACK_INCLUDE_TOPLEVEL_DIRECTORY 1)
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "MMG: 2d, surface and 3d remeshers")
SET(CPACK_PACKAGE_VENDOR "Cécile Dobrzynski, Pascal Frey, Charles Dapogny,"
  " Algiane Froehly")
SET(CPACK_PACKAGE_CONTACT "contact@mmgtools.org")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION ${CMAKE_RELEASE_VERSION})
SET(CPACK_PACKAGE_VERSION_MAJOR ${CMAKE_RELEASE_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${CMAKE_RELEASE_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${CMAKE_RELEASE_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME}-${CPACK_PACKAGE_VERSION}")
SET(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_INSTALL_DIRECTORY}-${CMAKE_SYSTEM})
SET(CPACK_PACKAGE_EXECUTABLES "${PROJECT_NAME}")
SET(CPACK_OUTPUT_FILE_PREFIX "../archives")

#SET(CPACK_SOURCE_IGNORE_FILES ";src/common/mmgtypes.h$" ${CPACK_SOURCE_IGNORE_FILE})
#SET(CPACK_SOURCE_IGNORE_FILES ";.\\\\\\\\.#;/\\\\\\\\.git;.~;.#;.-e;.DS_Store"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";lib/;bin;archives;def"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";.dir-locals.el;.tm_properties"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/mmg3d/example0/example0_a/mesh.o.*"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/mmg3d/example0/example0_b/mesh.o.*"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/mmg3d/example1/result0.*"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/mmg3d/example2/result1.*"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";libexamples/mmg3d/example3/SphereIso.*"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";CMakeFiles" ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";cmake/supp/libSystem.supp"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";CTestConfig.cmake" ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";cmake/testing/testing.cmake"
#  ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";.#$"      ${CPACK_SOURCE_IGNORE_FILES})
#SET(CPACK_SOURCE_IGNORE_FILES ";build*" ${CPACK_SOURCE_IGNORE_FILES}  )
#SET(CPACK_SOURCE_IGNORE_FILES ";doc/doxygen/.*/html;doc/doxygen/.*/Doxyfile$"
#  ${CPACK_SOURCE_IGNORE_FILES}  )

#IF ( WIN32 AND (NOT UNIX) )
#  #  # Warning compilation under windows not tested
#  #  #SET(CPACK_PACKAGE_ICON "${CMAKE_SOURCE_DIR}/logo.pdf")
#  #  SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\${PROJECT_NAME}.exe")
#  #  SET(CPACK_NSIS_DISPLAY_NAME
#  #    "${CPACK_PACKAGE_INSTALL_DIRECTORY} mmg: surface and volume remeshers")
#  #  SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\www.mmgtools.org‎")
#  #  SET(CPACK_NSIS_URL_INFO_ABOUT
#  #    "http:\\\\\\\\www.math.u-bordeaux1.fr/~dobrzyns/.../mmg3d.php‎")
#  #  SET(CPACK_NSIS_CONTACT ${CPACK_PACKAGE_CONTACT})
#  #  SET(CPACK_NSIS_MODIFY_PATH ON) # Choice to add prog dir to the system PATH
#  #  SET(CPACK_PACKAGING_INSTALL_PREFIX "${PROJECT_NAME}" CACHE PATH
#  #    "Where you want to install your package")
#  #  MARK_AS_ADVANCED(CPACK_PACKAGING_INSTALL_PREFIX)
#  #ELSEIF (APPLE )
#  #  SET(CPACK_BINARY_OSX11 "ON")
#  #ELSEIF ( UNIX )
#ENDIF (  WIN32 AND (NOT UNIX) )

INCLUDE(CPack)
