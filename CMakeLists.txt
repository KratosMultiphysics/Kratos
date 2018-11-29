project (KratosMultiphysics)
cmake_minimum_required (VERSION 2.8.6)

# Do not recursively expand variables
if(POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW)
endif(POLICY CMP0054)

# Uses INTERFACE_LINK_LIBRARIES instead of LINK_INTERFACE_LIBRARIES
if(POLICY CMP0022)
    cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)

# Set here the version number **** only update upon tagging a release!
set (KratosMultiphysics_MAJOR_VERSION 6)
set (KratosMultiphysics_MINOR_VERSION 1)
set (KratosMultiphysics_PATCH_VERSION Dev)

# If KRATOS_SOURCE_DIR is not defined use the CMAKE_SOURCE_DIR
if(NOT DEFINED KRATOS_SOURCE_DIR)
  message("-- No KRATOS_SOURCE_DIR is defined, using: ${CMAKE_SOURCE_DIR}")
  set (KRATOS_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
endif(NOT DEFINED KRATOS_SOURCE_DIR)

# Set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${KRATOS_SOURCE_DIR}/cmake_modules")

# Include cmake modules
include(DownloadLib)

# Define custom compiler build types
set( CMAKE_CONFIGURATION_TYPES Release RelWithDebInfo Debug FullDebug Custom)
set( BASIC_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
#set( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" )
#set( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
set( CMAKE_CXX_FLAGS_DEBUG "${BASIC_DEBUG_FLAGS} -DNDEBUG -DKRATOS_DEBUG" )
set( CMAKE_CXX_FLAGS_FULLDEBUG "${BASIC_DEBUG_FLAGS} -DKRATOS_DEBUG" )
#set( CMAKE_CXX_FLAGS_CUSTOM "${CMAKE_CXX_FLAGS_CUSTOM}" )

# Define internal CMake flags needed
set( CMAKE_C_FLAGS_FULLDEBUG "${CMAKE_C_FLAGS_DEBUG}" )
set( CMAKE_SHARED_LINKER_FLAGS_FULLDEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
set( CMAKE_EXE_LINKER_FLAGS_FULLDEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" )
set( CMAKE_C_FLAGS_CUSTOM "" )
set( CMAKE_CXX_FLAGS_CUSTOM "" )
set( CMAKE_SHARED_LINKER_FLAGS_CUSTOM "" )
set( CMAKE_EXE_LINKER_FLAGS_CUSTOM "" )
set( CMAKE_MODULE_LINKER_FLAGS_FULLDEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
set( CMAKE_MODULE_LINKER_FLAGS_CUSTOM "${CMAKE_MODULE_LINKER_FLAGS}")

# Search Git and the SHA1 Associated with the commit in the HEAD
find_package(Git)
if(GIT_FOUND)
    execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
        OUTPUT_VARIABLE KratosMultiphysics_SHA1_NUMBER
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else(GIT_FOUND)
    message("Git was not found on your system. SHA1 number will be set to 0.")
    set (KratosMultiphysics_SHA1_NUMBER 0)
endif(GIT_FOUND)

# Enable COTIRE
if(USE_COTIRE MATCHES ON)
    set(CMAKE_SKIP_RPATH ON)
endif(USE_COTIRE MATCHES ON)

# If build mode is not defined, assume Release
if(DEFINED ENV{KRATOS_BUILD_TYPE})
    list(FIND CMAKE_CONFIGURATION_TYPES $ENV{KRATOS_BUILD_TYPE} FOUND_MODE)
    if(FOUND_MODE EQUAL -1)
        message("Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}', using 'Release'.")
        set(CMAKE_BUILD_TYPE Release)
    else(FOUND_MODE EQUAL -1)
        set(CMAKE_BUILD_TYPE $ENV{KRATOS_BUILD_TYPE})
    endif(FOUND_MODE EQUAL -1)
else(DEFINED ENV{KRATOS_BUILD_TYPE})
    message("Undefined CMAKE_BUILD_TYPE, building in 'Release' mode.")
    set(CMAKE_BUILD_TYPE Release)
endif(DEFINED ENV{KRATOS_BUILD_TYPE})

# Configure files depending on the build type
set(KratosMultiphysics_BUILD_TYPE ${CMAKE_BUILD_TYPE})
if(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(KratosMultiphysics_BUILD_SUFFIX "")
else(${CMAKE_BUILD_TYPE} MATCHES "Release")
    set(KratosMultiphysics_BUILD_SUFFIX "-${CMAKE_BUILD_TYPE}")
endif(${CMAKE_BUILD_TYPE} MATCHES "Release")

configure_file(
    "${PROJECT_SOURCE_DIR}/kratos_version.h.in"
    "${PROJECT_SOURCE_DIR}/kratos/includes/kratos_version.h"
)

# Echo build type
message("-- The build type selected is: ${CMAKE_BUILD_TYPE}")

# Generate the export header for win.
include (GenerateExportHeader)

# Set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -Wall -std=c++11")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.0)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wsuggest-override")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -Wall")
    set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
endif(${CMAKE_COMPILER_IS_GNUCXX})

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -Wall -Wno-unknown-pragmas")
        set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        # pybind requires min. version 17, but we need at least 18:
            message( "DEPRECATED: detected compiler as Intel " ${CMAKE_CXX_COMPILER_VERSION} )
            message( FATAL_ERROR "Please use Version 18 or greater")
        endif()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops -lpthread -wd654 -wd10010")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -lpthread -wd654 -wd10010")
        set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
else(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops")
        set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Print the message with the custom compiler flags
message("-- Additional default options were set for ${CMAKE_CXX_COMPILER_ID} compiler:")
message("\tCMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message("\tCMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Old non-compatible versions of VS
if(${MSVC80})
    message( "DEPRECATED: detected compiler as MSVC80")
    message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC80})

if(${MSVC90})
    message( "DEPRECATED: detected compiler as MSVC90")
    message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC90})

if(${MSVC10})
    message( "DEPRECATED: detected compiler as MSVC10")
    message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC10})

if(${MSVC11})
    message( "DEPRECATED: detected compiler as MSVC11")
    message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC11})

if(${MSVC12})
    message( "DEPRECATED: detected compiler as MSVC12")
    message( FATAL_ERROR "Please use VisualStudio 2015 or greater")
endif(${MSVC12})

# Common flags for MSVC
if(${MSVC})
    message("Detected compiler as MSVC")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
    string( REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS} )
    string( REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
endif(${MSVC})

# Specific flags for different versions of MSVC
if(${MSVC14})
    # message("Adding additional flags for MSVC14")
    # Nothing needed right now
endif(${MSVC14})

# If no test policy enable by default
if(NOT DEFINED KRATOS_BUILD_TESTING)
    set (KRATOS_BUILD_TESTING ON)
endif(NOT DEFINED KRATOS_BUILD_TESTING)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    message("-- Kratos Tests are: Enabled")
else(${KRATOS_BUILD_TESTING} MATCHES ON)
    message("-- Kratos Tests are: Disabled")
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

# Pybind11
include(pybind11Tools)

if(USE_COTIRE MATCHES ON)
    include(cotire)
endif(USE_COTIRE MATCHES ON)

# Reset pybind11 config and remove -LTO on FullDebug to speedup linking time,
if(CMAKE_BUILD_TYPE MATCHES FullDebug)
    set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
    set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")
endif(CMAKE_BUILD_TYPE MATCHES FullDebug)

# Set installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Setting the cache path prevents it to change in case someone 'make' after modifying this file and not reconfiguring
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin/${CMAKE_SYSTEM_NAME}/${CMAKE_BUILD_TYPE}" CACHE PATH "Default Install path" FORCE)
    message("-- Standard install dir ${CMAKE_INSTALL_PREFIX}")
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message("-- User defined install dir ${CMAKE_INSTALL_PREFIX}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

# *****************************
# Find all dependencies
# *****************************

# Find and include OpenMP if possible
if(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
    message("****************************************************************************************")
    message("** OpenMP calculation was excluced artificially. Remove the var                       **")
    message("** KRATOS_EXCLUDE_OPENMP from the configuration script or set it to off to include it **")
    message("****************************************************************************************")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -diag-disable 3180 ") # Disabling warning abt unrecognized omp pragmas
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
else(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
    message("-- OpenMP included")
    find_package(OpenMP)
    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif (OPENMP_FOUND)
endif(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)

# Finding and including BOOST library (version should not matter anymore)
find_package(Boost)

if(Boost_FOUND)
  message("-- Boost found in ${BOOST_INCLUDE_DIRS}")
else(Boost_FOUND)
    DownloadLib("Boost" "https://dl.bintray.com/boostorg/release/1.67.0/source/boost_1_67_0.zip")
    set(BOOST_ROOT "${KRATOS_SOURCE_DIR}/external_libraries/boost_1_67_0")
    find_package(Boost)
endif(Boost_FOUND)

# If not found, download:

set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_REALPATH ON)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})

message("-- Boost Include: ${Boost_INCLUDE_DIRS}")
message("-- Boost Linkdir: ${Boost_LIBRARY_DIRS}")

# Find CUDA
if(${USE_CUDA} MATCHES ON)
    find_package(CUDA QUIET REQUIRED)
endif(${USE_CUDA} MATCHES ON)

# Follow the symbolic links for the python lib only if needed
if(${INSTALL_RUNKRATOS} MATCHES ON )
    foreach( library ${PYTHON_LIBRARIES})
        get_filename_component(python_lib ${library} REALPATH)
        set(PYTHON_LIBRARIES_REALPATH ${PYTHON_LIBRARIES_REALPATH} ${python_lib})
    endforeach(library)

    set(PYTHON_LIBRARIES ${PYTHON_LIBRARIES_REALPATH} CACHE STRING "exact file of the PYTHON library" FORCE)
endif(${INSTALL_RUNKRATOS} MATCHES ON )

# Enable BLAS and LAPACK
if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )
    # If prescribed set the fortran mangling
    if(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
        add_definitions(${KRATOS_SUPERLU_FORTRAN_MANGLING})
    else(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
        add_definitions(-DADD_)
    endif(DEFINED ${KRATOS_SUPERLU_FORTRAN_MANGLING})
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )

# Enable TETGEN
if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    find_package(Tetgen REQUIRED)
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

# Enable MPI
if(${METIS_APPLICATION} MATCHES ON )
    set(MPI_NEEDED ON)
endif(${METIS_APPLICATION} MATCHES ON )
if(${TRILINOS_APPLICATION} MATCHES ON )
    set(MPI_NEEDED ON)
endif(${TRILINOS_APPLICATION} MATCHES ON )
if(${MPI_SEARCH_APPLICATION} MATCHES ON )
    set(MPI_NEEDED ON)
endif(${MPI_SEARCH_APPLICATION} MATCHES ON )

# Find MPI
if(${MPI_NEEDED} MATCHES ON )
    find_package(MPI REQUIRED)
    add_definitions( -DKRATOS_USING_MPI )
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

    # Note: that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(${MPI_NEEDED} MATCHES ON )

# Find ZLIB (for gidpost)
find_package(ZLIB)

if( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
else( ZLIB_FOUND )
    message( "-- Preparing local ZLIB compilation.")
    include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/zlib )
    # Compile our own
    add_subdirectory(external_libraries/zlib)
    # we have to add this for zconf
    include_directories( ${CMAKE_BINARY_DIR}/external_libraries/zlib )
    set( ZLIB_LIBRARIES zlib )
endif( ZLIB_FOUND )

# including the AMatrix
if( DEFINED ENV{AMATRIX_DIR} )
    message("-- AMatrix found : ${AMATRIX_DIR}")
    include_directories ( "${AMATRIX_DIR}/include" )
    add_definitions(-DKRATOS_USE_AMATRIX)
else()
    message( "-- AMatrix not found. Please define the AMATRIX_DIR pointing to the AMatrix root directory.")
endif()

######################################################################################

# Include dir for external libraries
include_directories( ${KRATOS_SOURCE_DIR}/external_libraries )

# defines needed
add_definitions( -DKRATOS_PYTHON )

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# External librarie included subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(external_libraries/triangle)

add_subdirectory(kratos)

# Applications
# Add kratos applications by directory (Kratos apps, located in ./applications)
message("-- This kratos applications will be compiled:")
if(DEFINED ENV{KRATOS_INT_APP})
    foreach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
        message("\t${APPLICATION_PATH}")
    endforeach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
endif(DEFINED ENV{KRATOS_INT_APP})

# Add external applications by directory (Thrid party apps, located elsewhere)
message("-- This external applications will be compiled:")
if(DEFINED ENV{KRATOS_EXT_APP})
    foreach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        message("\t${APPLICATION_NAME}")
    endforeach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
endif(DEFINED ENV{KRATOS_EXT_APP})

# Configure kratos applications
message("-- Configuring kratos applications:")
if(DEFINED ENV{KRATOS_INT_APP})
    foreach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
        add_subdirectory(applications/${APPLICATION_PATH})
    endforeach(APPLICATION_PATH $ENV{KRATOS_INT_APP})
endif(DEFINED ENV{KRATOS_INT_APP})

# Configure kratos applications
message("-- Configuring external applications:")
if(DEFINED ENV{KRATOS_EXT_APP})
    foreach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        add_subdirectory(${APPLICATION_PATH} ${CMAKE_CURRENT_BINARY_DIR}/applications/${APPLICATION_NAME})
    endforeach(APPLICATION_PATH $ENV{KRATOS_EXT_APP})
endif(DEFINED ENV{KRATOS_EXT_APP})

# Embedded python - compiles "krun and packages the python lib needed"
if(${INSTALL_RUNKRATOS} MATCHES ON )
    add_subdirectory(runkratos)
endif(${INSTALL_RUNKRATOS} MATCHES ON )

# Mpi
if(${MPI_NEEDED} MATCHES ON )
    add_subdirectory(external_libraries/mpi_python)
endif(${MPI_NEEDED} MATCHES ON )

# Install Files
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/__init__.py" DESTINATION KratosMultiphysics )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/kratos_utilities.py" DESTINATION KratosMultiphysics )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" DESTINATION KratosMultiphysics )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/application_importer.py" DESTINATION KratosMultiphysics )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/kratos_unittest.py" DESTINATION KratosMultiphysics RENAME KratosUnittest.py )
install(FILES "${KRATOS_SOURCE_DIR}/kratos/python_interface/KratosLoader.py" DESTINATION KratosMultiphysics RENAME KratosLoader.py )

# Install the libraries in the libs folder
install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)

# Kratos Testing. Install everything except sources
if(${INSTALL_TESTING_FILES} MATCHES ON )
    install(DIRECTORY ${KRATOS_SOURCE_DIR}/kratos/tests/ DESTINATION ${CMAKE_INSTALL_PREFIX}/kratos/tests
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c"   EXCLUDE
        PATTERN "*.h"   EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)