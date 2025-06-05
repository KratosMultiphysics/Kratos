cmake_minimum_required (VERSION 3.15.0)

project (KratosMultiphysics)
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# Defining the base source folder
set(KRATOS_BASE_FOLDER "${CMAKE_CURRENT_SOURCE_DIR}")

# Setting some policies
# No recursive dereferencing
if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif(POLICY CMP0054)

# Uses INTERFACE_LINK_LIBRARIES instead of LINK_INTERFACE_LIBRARIES
if(POLICY CMP0022)
  cmake_policy(SET CMP0022 NEW)
endif(POLICY CMP0022)

# Find package based on <packagename>_ROOT variable (both policies are needed)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif(POLICY CMP0074)

if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif(POLICY CMP0144)

# Lists
if(POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif(POLICY CMP0057)

# Timestamp for extracted files
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif(POLICY CMP0135)

# Path Normalization
if(POLICY CMP0177)
  cmake_policy(SET CMP0177 NEW)
endif(POLICY CMP0177)

# Boost config
if(POLICY CMP0167)
  cmake_policy(SET CMP0167 OLD)
endif(POLICY CMP0167)

# Set here the version number **** only update upon tagging a release!
set (KratosMultiphysics_MAJOR_VERSION 10)
set (KratosMultiphysics_MINOR_VERSION 2)
set (KratosMultiphysics_PATCH_VERSION 1)

# If KRATOS_SOURCE_DIR is not defined use the CMAKE_SOURCE_DIR
if(NOT DEFINED KRATOS_SOURCE_DIR)
  message("-- No KRATOS_SOURCE_DIR is defined, using: ${CMAKE_SOURCE_DIR}")
  set (KRATOS_SOURCE_DIR "${CMAKE_SOURCE_DIR}")
endif(NOT DEFINED KRATOS_SOURCE_DIR)

if(NOT DEFINED KRATOS_BINARY_DIR)
  message("-- No KRATOS_BINARY_DIR is defined, using: ${CMAKE_BINARY_DIR}")
  set(KRATOS_BINARY_DIR ${CMAKE_BINARY_DIR})
endif(NOT DEFINED KRATOS_BINARY_DIR)

option(REMOVE_INSTALL_DIRECTORIES "Remove the install directories to achieve a clean compilation" ON)

# Link on install
if(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
  message("-- No KRATOS_INSTALL_PYTHON_USING_LINKS is defined, setting to OFF")
  SET(INSTALL_PYTHON_USING_LINKS OFF)
else(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
  SET(INSTALL_PYTHON_USING_LINKS $ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})
endif(NOT DEFINED ENV{KRATOS_INSTALL_PYTHON_USING_LINKS})

# Setting the libs folder for the shared objects built in kratos
if (APPLE)
  cmake_policy(SET CMP0042 NEW)
  set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/libs")
else()
  set(CMAKE_INSTALL_RPATH "$ORIGIN/../libs")
endif()

# If no test policy enable by default
option(KRATOS_BUILD_TESTING "KRATOS_BUILD_TESTING defines if the C++ tests are compiled. These increase compilation time, but if not set only python tests can be executed. Default setting is ON" ON)

# If no benchmark policy enable by default
option(KRATOS_BUILD_BENCHMARK "KRATOS_BUILD_BENCHMARK defines if the C++ benchmarks (Google benchmark) are compiled. These increase compilation time. Default setting is OFF" OFF)

# If logger coloring policy enable by default
option(KRATOS_COLORED_OUTPUT "KRATOS_COLORED_OUTPUT defines if the logger output it is colored. Default setting is OFF" OFF)

# If no pch policy disable by default
option(KRATOS_USE_PCH "KRATOS_USE_PCH defines if pch will be used during the compilation. This may will decrease compilation for clean build or if core components are touched. Default setting is OFF" OFF)

# If not defined mute the very verbose installation messages
if(NOT DEFINED CMAKE_INSTALL_MESSAGE)
  message("-- No CMAKE_INSTALL_MESSAGE is defined. Setting to NEVER")
  SET (CMAKE_INSTALL_MESSAGE NEVER)
endif(NOT DEFINED CMAKE_INSTALL_MESSAGE)

# Set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${KRATOS_SOURCE_DIR}/cmake_modules")

# Generate and copy configuration for language servers
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Include cmake modules
include(DownloadLib)
include(ExternalProject)

# Define custom compiler build types
SET( CMAKE_CONFIGURATION_TYPES Release RelWithDebInfo Debug FullDebug Custom)
SET( BASIC_DEBUG_FLAGS "${CMAKE_CXX_FLAGS_DEBUG}")
#SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" )
#SET( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" )
SET( CMAKE_CXX_FLAGS_DEBUG "${BASIC_DEBUG_FLAGS} -DNDEBUG -DKRATOS_DEBUG" )
SET( CMAKE_CXX_FLAGS_FULLDEBUG "${BASIC_DEBUG_FLAGS} -D_DEBUG -DKRATOS_DEBUG" )
SET( CMAKE_CXX_FLAGS_CUSTOM "${CMAKE_CXX_FLAGS_CUSTOM}" )
set( CMAKE_MAP_IMPORTED_CONFIG_FULLDEBUG Debug)

# Define internal CMake flags needed
SET( CMAKE_C_FLAGS_FULLDEBUG "${CMAKE_C_FLAGS_DEBUG}" )
SET( CMAKE_SHARED_LINKER_FLAGS_FULLDEBUG "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_EXE_LINKER_FLAGS_FULLDEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}" )
SET( CMAKE_C_FLAGS_CUSTOM "" )
SET( CMAKE_CXX_FLAGS_CUSTOM "" )
SET( CMAKE_SHARED_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_EXE_LINKER_FLAGS_CUSTOM "" )
SET( CMAKE_MODULE_LINKER_FLAGS_FULLDEBUG "${CMAKE_MODULE_LINKER_FLAGS}")
SET( CMAKE_MODULE_LINKER_FLAGS_CUSTOM "${CMAKE_MODULE_LINKER_FLAGS}")

# If build mode is not defined, assume Release
if(NOT CMAKE_BUILD_TYPE)
  if(DEFINED ENV{KRATOS_BUILD_TYPE})
      list(FIND CMAKE_CONFIGURATION_TYPES $ENV{KRATOS_BUILD_TYPE} FOUND_MODE)
      if(FOUND_MODE EQUAL -1)
          message("Unknown CMAKE_BUILD_TYPE '${CMAKE_BUILD_TYPE}', using 'Release'.")
          set(CMAKE_BUILD_TYPE Release)
      else(FOUND_MODE EQUAL -1)
          set(CMAKE_BUILD_TYPE $ENV{KRATOS_BUILD_TYPE})
      endif(FOUND_MODE EQUAL -1)
  else(DEFINED ENV{KRATOS_BUILD_TYPE})
      message("--No CMAKE_BUILD_TYPE is defined, building in 'Release' mode.")
      set(CMAKE_BUILD_TYPE Release)
  endif(DEFINED ENV{KRATOS_BUILD_TYPE})
endif(NOT CMAKE_BUILD_TYPE)
set (KratosMultiphysics_BUILD_TYPE ${CMAKE_BUILD_TYPE})

# Unitary builds with version is lower than 3.16 deactivate
IF(CMAKE_VERSION VERSION_LESS "3.16.0")
    IF(CMAKE_UNITY_BUILD MATCHES ON)
        SET(CMAKE_UNITY_BUILD OFF)
    ENDIF(CMAKE_UNITY_BUILD MATCHES ON)
ELSE(CMAKE_VERSION VERSION_LESS "3.16.0")
    IF(USE_COTIRE MATCHES ON)
        # Legacy Unitary build
        MESSAGE(STATUS "Using USE_COTIRE is deprecated and was replaced by CMAKE_UNITY_BUILD.")
        SET(CMAKE_UNITY_BUILD ON)
    ENDIF(USE_COTIRE MATCHES ON)
    # Some checks
    IF(CMAKE_UNITY_BUILD MATCHES ON)
        IF(NOT DEFINED KRATOS_UNITY_BUILD_BATCH_SIZE)
            SET(KRATOS_UNITY_BUILD_BATCH_SIZE 32)
        ENDIF(NOT DEFINED KRATOS_UNITY_BUILD_BATCH_SIZE)
        IF(${MSVC})
          SET(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY True)
        ENDIF(${MSVC})
    ENDIF(CMAKE_UNITY_BUILD MATCHES ON)
ENDIF(CMAKE_VERSION VERSION_LESS "3.16.0")

# Get subversion data. This is done automagically by the cmakes
include (GenerateExportHeader)

# Search the SHA1 Associated with the commit in the HEAD
find_package(Git)
if(GIT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    OUTPUT_VARIABLE KratosMultiphysics_SHA1_NUMBER
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE SHA1_NOT_FOUND
  )
  if(SHA1_NOT_FOUND)
    message("Git did not find the SHA1 number. It will be set to 0.")
    set (KratosMultiphysics_SHA1_NUMBER 0)
  endif(SHA1_NOT_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} branch --show-current
    OUTPUT_VARIABLE KratosMultiphysics_BRANCH_NAME
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE BRANCH_NAME_NOT_FOUND
  )
  if(BRANCH_NAME_NOT_FOUND)
    message("Git did not find the branch name. It will be set to empty.")
    set (KratosMultiphysics_BRANCH_NAME "")
  endif(BRANCH_NAME_NOT_FOUND)
else(GIT_FOUND)
  message(STATUS "Git was not found on your system. SHA1 number will be set to 0. Branch name to empty.")
  set (KratosMultiphysics_SHA1_NUMBER 0)
  set (KratosMultiphysics_BRANCH_NAME "")
endif(GIT_FOUND)

# Configure files depending on the build type
if( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "")
else( ${CMAKE_BUILD_TYPE} MATCHES "Release" )
  set (KratosMultiphysics_BUILD_SUFFIX "-${CMAKE_BUILD_TYPE}" )
endif( ${CMAKE_BUILD_TYPE} MATCHES "Release" )

# Set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -Wall -Wl,--no-as-needed -ldl -Wsuggest-override")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 8.0)
    message(STATUS "DEPRECATED: detected compiler as GCC " ${CMAKE_CXX_COMPILER_VERSION} )
    message(FATAL_ERROR "Please use Version 8 or greater")
  endif()
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -Wall -Wl,--no-as-needed -ldl")
  set(CMAKE_SHARED_LINKER_FLAGS " -Wl,--no-undefined")
  message(STATUS "Additional default options were set for GCC")
endif(${CMAKE_COMPILER_IS_GNUCXX})

# Check for correct clang version
if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
    message( FATAL_ERROR "Please use Version 12.0 or greater")
endif()

message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
      if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas  ")
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -funroll-loops -Wall -Wno-unknown-pragmas  ")
        set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
        message(STATUS "Additional default options were set for LLVM/Clang compiler")
      endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -funroll-loops -Wall -Wno-unused-local-typedef -Wno-unknown-pragmas ")
      set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops -Wall -Wno-unknown-pragmas ")
      set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
      message(STATUS "Additional default options were set for LLVM/Clang compiler")
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    if(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
      if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 18.0)
        # pybind requires min. version 17, but we need at least 18:
        message(STATUS "DEPRECATED: detected compiler as Intel " ${CMAKE_CXX_COMPILER_VERSION} )
        message(FATAL_ERROR "Please use Version 18 or greater")
      endif()
      if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        # "IntelLLVM": Intel's NextGen LLVM compiler (icx)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -funroll-loops -diag-disable 654 -diag-disable 10010 -diag-disable 1011")
      else (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
        # "Intel": Intel's Classic compiler (icc)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -funroll-loops -lpthread -wd654 -wd10010 -wd1011")
      endif()
      set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  ")
      set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
      message(STATUS "Additional default options were set for Intel Compiler")
    endif(${CMAKE_CXX_COMPILER_ID} MATCHES Intel)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
else(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  if(${CMAKE_C_COMPILER} MATCHES "icc.*$")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC  -funroll-loops  ")
    set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  ")
    set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")
    message(STATUS "Additional default options were set for Intel Compiler")
  endif(${CMAKE_C_COMPILER} MATCHES "icc.*$")
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# After all the compiler check we can print the final flags
message(STATUS "CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
message(STATUS "CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")

# Old non-compatible versions of VS
if(${MSVC})
  if(${MSVC_TOOLSET_VERSION} LESS 141)
    message(STATUS "DEPRECATED: detected compiler as ${MSVC_TOOLSET_VERSION}")
    message(FATAL_ERROR "Please use VisualStudio Toolset 141 (VisualStudio 2017) or greater")
  endif(${MSVC_TOOLSET_VERSION} LESS 141)
endif(${MSVC})

# Common flags for MSVC
if(${MSVC})
  message(STATUS "Detected compiler as MSVC")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W1 /bigobj /EHsc -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W1 /bigobj /EHsc /Zc:__cplusplus -DBOOST_ALL_NO_LIB -D_SCL_SECURE_NO_WARNINGS")
  string( REPLACE "/W3" "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS} )
  string( REPLACE "/W3" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} )
  # CMake does not select the correct runtime library for FullDebug now, no idea why
  if( ${CMAKE_BUILD_TYPE} MATCHES "FullDebug" )
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDebugDLL" CACHE STRING "" FORCE)
  endif( ${CMAKE_BUILD_TYPE} MATCHES "FullDebug" )
endif(${MSVC})

# Specific flags for different versions of MSVC
if(${MSVC14})
  # message("Adding additional flags for MSVC14")
  # Nothing needed right now
endif(${MSVC14})

# If install prefix is set, install all the files.
SET(INSTALL_TESTING_FILES ON)
SET(INSTALL_PYTHON_FILES ON)        # To be removed when all applications are ported

#include the file with custom installation properties
include(install_function)
include(KratosDependencies)
include(KratosGTest)
include(KratosAddSources)
include(KratosGBenchmark)
include(FetchContent)

# Logger configuration
if(KRATOS_COLORED_OUTPUT MATCHES ON)
  add_definitions(-DKRATOS_COLORED_OUTPUT)
endif(KRATOS_COLORED_OUTPUT MATCHES ON)

# Testing
if(KRATOS_BUILD_TESTING MATCHES ON)
  # Add the definitions if required
  ADD_DEFINITIONS(-DKRATOS_BUILD_TESTING)

  # retrieve a copy of the current directory's `COMPILE_OPTIONS`
  get_directory_property(kratos_root_compile_options COMPILE_OPTIONS)

  # Disable warnings (needed by Intel Compiler Legacy)
  add_compile_options(-w)

  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
  # GTest doc is wrong, this option should be used
  set(BUILD_SHARED_LIBS ON)
  set(CMAKE_INSTALL_LIBDIR "libs") # In Linux default dir is lib
  set(CMAKE_INSTALL_BINDIR "libs") # In Win default dir is bin
  FetchContent_MakeAvailable(googletest)
  enable_testing()

  # restore the current directory's old `COMPILE_OPTIONS`
  set_directory_properties(PROPERTIES COMPILE_OPTIONS "${kratos_root_compile_options}")
endif(KRATOS_BUILD_TESTING MATCHES ON) 

# Benchmarking
if(KRATOS_BUILD_BENCHMARK MATCHES ON)
  # Add the definitions if required
  ADD_DEFINITIONS(-DKRATOS_BUILD_BENCHMARKING)

  # Retrieve a copy of the current directory's `COMPILE_OPTIONS`
  get_directory_property(kratos_root_compile_options COMPILE_OPTIONS)

  # Disable warnings (needed by centos. We should all love centos, it clearly needs some affection)
  add_compile_options(-w)

  FetchContent_Declare(
    googlebenchmark
    URL https://github.com/google/benchmark/archive/v1.9.0.zip
  )
  # For Windows: Prevent overriding the parent project's compiler/linker settings
  set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)
  set(CMAKE_INSTALL_LIBDIR "libs") # In Linux default dir is lib
  set(CMAKE_INSTALL_BINDIR "libs") # In Win default dir is bin
  FetchContent_MakeAvailable(googlebenchmark)

  # Restore the current directory's old `COMPILE_OPTIONS`
  set_directory_properties(PROPERTIES COMPILE_OPTIONS "${kratos_root_compile_options}")
endif(KRATOS_BUILD_BENCHMARK MATCHES ON)

################### PYBIND11

# Try to use python executable from env variable
if(DEFINED ENV{PYTHON_EXECUTABLE})
  set(PYTHON_EXECUTABLE $ENV{PYTHON_EXECUTABLE})
endif(DEFINED ENV{PYTHON_EXECUTABLE})

include(pybind11Tools)

# Reset pybind11 config and remove -LTO since it gives multiple problems.
if(NOT DEFINED KRATOS_ENABLE_LTO)
    message(STATUS "LTO is Disabled")
    set(PYBIND11_LTO_CXX_FLAGS "" CACHE INTERNAL "")
    set(PYBIND11_LTO_LINKER_FLAGS "" CACHE INTERNAL "")
else(NOT DEFINED KRATOS_ENABLE_LTO)
    message(STATUS "LTO is Enabled")
    if(${CMAKE_COMPILER_IS_GNUCXX})
        set(CMAKE_AR "gcc-ar")
        set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> qcs <TARGET> <LINK_FLAGS> <OBJECTS>")
        set(CMAKE_C_ARCHIVE_FINISH true)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto=jobserver")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=jobserver")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto=jobserver")
    elseif(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -flto=thin")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=thin")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -flto=thin")
    endif()
endif(NOT DEFINED KRATOS_ENABLE_LTO)

message(STATUS "AR VERSION: ${CMAKE_AR}")

# check version of Python, needs to be done after including pybind
if(${PYTHON_VERSION_MAJOR} LESS 3 OR (${PYTHON_VERSION_MAJOR} EQUAL 3 AND ${PYTHON_VERSION_MINOR} LESS 8))
    message( FATAL_ERROR "Kratos only supports Python version 3.8 and above")
endif()
set(PYTHON_INTERFACE_VERSION "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

message(STATUS "Python version used for the interface will be ${PYTHON_INTERFACE_VERSION}")

# Set installation directory. TODO: Delete this and use CMAKE_INSTALL_PREFIX
if(DEFINED KRATOS_INSTALL_PREFIX)
    set(CMAKE_INSTALL_PREFIX ${KRATOS_INSTALL_PREFIX} )
endif(DEFINED KRATOS_INSTALL_PREFIX)

# Set installation directory
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Setting the cache path prevents it to change in case someone 'make' after modifying this file and not reconfiguring
    set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/bin/${CMAKE_BUILD_TYPE}" CACHE PATH "Default Install path" FORCE)
    message(STATUS "Standard install dir ${CMAKE_INSTALL_PREFIX}")
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    message(STATUS "User defined install dir ${CMAKE_INSTALL_PREFIX}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

## Find all dependencies

# Select shared memory parallelization
if(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)
    message(FATAL_ERROR "The option \"KRATOS_EXCLUDE_OPENMP\" was removed, use KRATOS_SHARED_MEMORY_PARALLELIZATION=\"None\" instead")
endif(${KRATOS_EXCLUDE_OPENMP} MATCHES ON)

# default to OpenMP
if(NOT DEFINED KRATOS_SHARED_MEMORY_PARALLELIZATION)
  message(STATUS "\"KRATOS_SHARED_MEMORY_PARALLELIZATION\" not defined, defaulting to \"OpenMP\"")
  SET(KRATOS_SHARED_MEMORY_PARALLELIZATION "OpenMP")
endif(NOT DEFINED KRATOS_SHARED_MEMORY_PARALLELIZATION)

if (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "OpenMP")
  find_package(OpenMP)
  if (OPENMP_FOUND)
    message(STATUS "Using OpenMP for shared memory parallelization")
    add_definitions( -DKRATOS_SMP_OPENMP )
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  else (OPENMP_FOUND)
    message(FATAL_ERROR "OpenMP could not be found!")
    # fallback solution => in future once better supported we can use the C++11 based parallelization instead
    message("OpenMP could not be found, disabling shared memory parallelization!")
    add_definitions( -DKRATOS_SMP_NONE )
  endif (OPENMP_FOUND)

elseif (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "C++11")
  message(STATUS "Using C++11 for shared memory parallelization")
  add_definitions( -DKRATOS_SMP_CXX11 )
  message(WARNING "Using C++11 for shared memory parallelization is highly experimental and not fully supported!")

elseif (KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "None")
  add_definitions( -DKRATOS_SMP_NONE )
  message(STATUS "Shared memory parallelization is disabled!")

else()
  message(FATAL_ERROR "Invalid option chosen for \"KRATOS_SHARED_MEMORY_PARALLELIZATION\"! Available options are: \"OpenMP\", \"C++11\", \"None\"")
endif ()

OPTION(KRATOS_NO_TRY_CATCH "Disabled the use of try/catch through the code" OFF)
if( KRATOS_NO_TRY_CATCH MATCHES ON )
  add_definitions( -DKRATOS_NO_TRY_CATCH )
endif( KRATOS_NO_TRY_CATCH MATCHES ON )

# Boost
include(KratosBoost)

# std::atomic_ref is part of C++20, hence using boost::atomic_ref
# boost::atomic_ref has full support from version 1.74
if((KRATOS_SHARED_MEMORY_PARALLELIZATION STREQUAL "C++11") AND (CMAKE_CXX_STANDARD LESS 20) AND (Boost_VERSION_STRING VERSION_LESS 1.74))
  message(FATAL_ERROR "Kratos requires at least boost version 1.74 when using C++11 for shared memory parallelization and compiling with C++17")
endif ()

##*****************************
# Finding cuda if needed
if(${USE_CUDA} MATCHES ON)
  find_package(CUDA QUIET REQUIRED)
endif(${USE_CUDA} MATCHES ON)

##*****************************
# Compiling the triangle library
if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    add_definitions(-DUSE_TRIANGLE_NONFREE_TPL)
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )

# Compiling the tetgen library
if(${USE_TETGEN_NONFREE_TPL} MATCHES ON)
    set(TETGEN_ZIP_FILE "tetgen")

    if(DEFINED USE_TETGEN_NONFREE_TPL_PATH)
        # Use TetGen from local source
        set(TETGEN_EXT_PATH ${CMAKE_CURRENT_BINARY_DIR}/external_libraries/tetgen)

        find_path(TETGEN_ROOT NAMES tetgen.h HINTS ${TETGEN_EXT_PATH}/*/)
        message("TETGEN_ROOT: ${TETGEN_ROOT}")
        if(NOT TETGEN_ROOT OR FORCE_TETGEN_NONFREE_TPL_URL)
            file(COPY ${USE_TETGEN_NONFREE_TPL_PATH} DESTINATION ${TETGEN_EXT_PATH})
            find_path(TETGEN_ROOT NAMES tetgen.h HINTS ${TETGEN_EXT_PATH}/*/)
            message("Copying tetgen dir to ${TETGEN_EXT_PATH}")
        endif(NOT TETGEN_ROOT OR FORCE_TETGEN_NONFREE_TPL_URL)

        # Prepare tetgen
        add_definitions(-DTETLIBRARY)
        add_definitions(${USE_TETGEN_NONFREE_TPL_FLAGS})

        set(Tetgen_FOUND True)

        # This is for backwards compatibility
        if(NOT DEFINED TETGEN_INCLUDE)
            set(TETGEN_INCLUDE ${TETGEN_ROOT})
        endif(NOT DEFINED TETGEN_INCLUDE)
        if(NOT DEFINED TETGEN_LIBRARIES)
            set(TETGEN_LIBRARIES tetgen)
        endif(NOT DEFINED TETGEN_LIBRARIES)
    elseif(DEFINED USE_TETGEN_NONFREE_TPL_URL)
        # Use TetGen from external source (download)
        set(TETGEN_EXT_PATH ${KRATOS_BASE_FOLDER}/external_libraries/tetgen)

        find_path(TETGEN_ROOT NAMES tetgen.h HINTS ${TETGEN_EXT_PATH}/*/)
        if(NOT TETGEN_ROOT OR FORCE_TETGEN_NONFREE_TPL_URL)
            DownloadLibAt(${TETGEN_ZIP_FILE} ${USE_TETGEN_NONFREE_TPL_URL} ${TETGEN_EXT_PATH})
            find_path(TETGEN_ROOT NAMES tetgen.h HINTS ${TETGEN_EXT_PATH}/*/)
        else(NOT TETGEN_ROOT OR FORCE_TETGEN_NONFREE_TPL_URL)
            # Touch the url variable so it does no trigger a unused warning.
            set(ignoreMe "${USE_TETGEN_NONFREE_TPL_URL}")
            message("-- Found tetgen at ${TETGEN_EXT_PATH}, skipping...")
        endif(NOT TETGEN_ROOT OR FORCE_TETGEN_NONFREE_TPL_URL)

        # Prepare tetgen
        add_definitions(-DTETLIBRARY)
        add_definitions(${USE_TETGEN_NONFREE_TPL_FLAGS})

        set(Tetgen_FOUND True)

        # This is for backwards compatibility
        if(NOT DEFINED TETGEN_INCLUDE)
            set(TETGEN_INCLUDE ${TETGEN_ROOT})
        endif(NOT DEFINED TETGEN_INCLUDE)
        if(NOT DEFINED TETGEN_LIBRARIES)
            set(TETGEN_LIBRARIES tetgen)
        endif(NOT DEFINED TETGEN_LIBRARIES)
    else(DEFINED USE_TETGEN_NONFREE_TPL_PATH)
        # Try to find tetgen installed in the system
        find_package(Tetgen REQUIRED)
    endif(DEFINED USE_TETGEN_NONFREE_TPL_PATH)

    if(tetgen_FOUND)
        message(FATAL_ERROR "TetGen could not be found. Please use -DUSE_TETGEN_NONFREE_TPL_PATH=[path] or -DUSE_TETGEN_NONFREE_TPL_URL=[url] to indicate the correct TetGen location.")
    endif(tetgen_FOUND)
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON)

# Setting MPI
if(${USE_MPI} MATCHES ON )
  find_package(MPI REQUIRED)
  # suppress annoying warnings from trilinos due to not having the same package name as the lib name.
  if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
  endif()
  find_package(TRILINOS QUIET)
  if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 0 CACHE INTERNAL "No dev warnings")
  endif()
  add_definitions( -DKRATOS_USING_MPI )
  set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
  set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

  # Note: that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
endif(${USE_MPI} MATCHES ON )

# If we do not locally compile ZLIB we will find ZLIB in the system, if we cannot find it we will enforce the local compilation
OPTION ( FORCE_LOCAL_ZLIB_COMPILATION "Force local compilation of zlib" OFF )
if( FORCE_LOCAL_ZLIB_COMPILATION MATCHES OFF )
    # Find zlib (for gidpost)
    find_package(ZLIB)
    if( ZLIB_FOUND )
        include_directories( ${ZLIB_INCLUDE_DIRS} )
    else( ZLIB_FOUND )
        set( FORCE_LOCAL_ZLIB_COMPILATION ON )
    endif( ZLIB_FOUND )
endif( FORCE_LOCAL_ZLIB_COMPILATION MATCHES OFF )

# We compile ZLIB locally
if( FORCE_LOCAL_ZLIB_COMPILATION MATCHES ON )
  message(STATUS "Preparing local ZLIB compilation.")
  include_directories( ${KRATOS_SOURCE_DIR}/external_libraries/zlib )
  # Compile our own
  add_subdirectory( ${KRATOS_SOURCE_DIR}/external_libraries/zlib )
  # We have to add this for zconf
  include_directories( ${KRATOS_BINARY_DIR}/external_libraries/zlib )
  set( ZLIB_LIBRARIES zlib )
endif( FORCE_LOCAL_ZLIB_COMPILATION MATCHES ON )

## Echo user options
message("\n-- CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}\n")

# Include dir for external libraries
include_directories( SYSTEM ${KRATOS_SOURCE_DIR}/external_libraries )

# defines needed
add_definitions( -DKRATOS_PYTHON )

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

set(${KRATOS_KERNEL} "")
set(${KRATOS_PYTHON_INTERFACE} "")

# Enable profiling if requested
if (${KRATOS_ENABLE_PROFILING} MATCHES ON)
    add_compile_definitions(KRATOS_ENABLE_PROFILING)
endif()

# Build exclusions and linking dirs.
# The intended usage of this flags is to allow a deferred/incremental compilation of core/applications.
# Please do not use this options for regular development and never use this with REMOVE_INSTALL_DIRECTORIES=ON
OPTION ( EXCLUDE_KRATOS_CORE "Force exclusion of the core" OFF )
OPTION ( EXCLUDE_AUTOMATIC_DEPENDENCIES "Force exclusion of automatically added app dependencies" OFF )

if(${REMOVE_INSTALL_DIRECTORIES} MATCHES ON AND (${EXCLUDE_KRATOS_CORE} MATCHES ON OR ${EXCLUDE_AUTOMATIC_DEPENDENCIES} MATCHES ON))
  message(FATAL_ERROR "REMOVE_INSTALL_DIRECTORIES cannot be enabled with EXCLUDE_KRATOS_CORE or EXCLUDE_AUTOMATIC_DEPENDENCIES")
endif()

if(${EXCLUDE_KRATOS_CORE} MATCHES ON OR ${EXCLUDE_AUTOMATIC_DEPENDENCIES} MATCHES ON)
  link_directories( ${CMAKE_INSTALL_PREFIX}/libs )
endif()

# Compiling the gismo library
add_subdirectory(${KRATOS_SOURCE_DIR}/external_libraries/gismo)
link_directories   (${GISMO_LIBRARY_DIR})


# Compile Kratos Core
if(${EXCLUDE_KRATOS_CORE} MATCHES OFF)
  add_subdirectory(external_libraries/gidpost)
  add_subdirectory(kratos)
endif()

# Configure kratos applications
message(STATUS "Configuring applications (ENV):")
set_property(GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES)
if(DEFINED ENV{KRATOS_APPLICATIONS})
    foreach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        list(APPEND LIST_OF_APPLICATIONS_TO_BE_COMPILED ${APPLICATION_NAME})
    endforeach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
endif(DEFINED ENV{KRATOS_APPLICATIONS})

if(DEFINED ENV{KRATOS_APPLICATIONS})
    foreach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
        get_filename_component(APPLICATION_NAME ${APPLICATION_PATH} NAME)
        if(NOT TARGET Kratos${APPLICATION_NAME})
            message("\tAdding application '${APPLICATION_PATH}'")
            add_subdirectory(${APPLICATION_PATH} ${CMAKE_CURRENT_BINARY_DIR}/applications/${APPLICATION_NAME})
        else(NOT TARGET Kratos${APPLICATION_NAME})
            message(WARNING "[Warning] Application '${APPLICATION_PATH}' was already added")
        endif(NOT TARGET Kratos${APPLICATION_NAME})
    endforeach(APPLICATION_PATH $ENV{KRATOS_APPLICATIONS})
endif(DEFINED ENV{KRATOS_APPLICATIONS})

# Create custom targets. Not the cleanest solution, but minimizes the changes
add_custom_target(KratosKernel DEPENDS ${KRATOS_KERNEL})
add_custom_target(KratosPythonInterface DEPENDS ${KRATOS_PYTHON_INTERFACE})

message("\n***********************************************************************\n")
message("        LIST OF APPLICATIONS THAT ARE GOING TO BE COMPILED:\n")
foreach(individual_app_name ${LIST_OF_APPLICATIONS_TO_BE_COMPILED})
  message("        ${individual_app_name}")
endforeach()
message("\n        APPLICATIONS ADDED TO COMPILATION THROUGH DEPENDENCIES:\n")
get_property(local_list GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES)
foreach(individual_app_name ${local_list})
    if(NOT ";${LIST_OF_APPLICATIONS_TO_BE_COMPILED};" MATCHES ${individual_app_name})
        message("        ${individual_app_name}")
    endif(NOT ";${LIST_OF_APPLICATIONS_TO_BE_COMPILED};" MATCHES ${individual_app_name})
endforeach()
message("\n***********************************************************************\n")

message("\nList of upcoming API-breaking and behavior changes:")
message("\t Geometry 'Dimension' removed: https://github.com/KratosMultiphysics/Kratos/pull/10977")
message("\t Geometry old projection methods deprecation: https://github.com/KratosMultiphysics/Kratos/pull/9024")
message("\t Geometry old closest point methods deprecation: https://github.com/KratosMultiphysics/Kratos/pull/9243")
message("\t Model 'CreateModelPart' behavior change (no error thrown if model part exists): https://github.com/KratosMultiphysics/Kratos/pull/9598")
message("\t Make GetIntegrationMethod method of Condition const: https://github.com/KratosMultiphysics/Kratos/pull/9769")
message("\t Behavior change of ModelPart.GetProperties (MeshIndex removed): https://github.com/KratosMultiphysics/Kratos/pull/9774")
message("\t Behavior change of Testing. Please ensure you use KRATOS_EXPECT for testing instead of KRATOS_CHECK")

# Compiling the clipper library
add_subdirectory(${KRATOS_SOURCE_DIR}/external_libraries/clipper)
# Compiling the triangle library
if(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )
    add_subdirectory(${KRATOS_SOURCE_DIR}/external_libraries/triangle)
endif(${USE_TRIANGLE_NONFREE_TPL} MATCHES ON )

# Compiling the tinyexpr library
add_subdirectory(${KRATOS_SOURCE_DIR}/external_libraries/tinyexpr)

# # Clean the Module and libs install directories
# if(${REMOVE_INSTALL_DIRECTORIES} MATCHES ON )
#     install(CODE "message(STATUS \"Deleting: ${CMAKE_INSTALL_PREFIX}/KratosMultiphysics\")")
#     install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/KratosMultiphysics\")")
#     install(CODE "message(STATUS \"Deleting: ${CMAKE_INSTALL_PREFIX}/libs\")")
#     install(CODE "file(REMOVE_RECURSE \"${CMAKE_INSTALL_PREFIX}/libs\")")
# endif(${REMOVE_INSTALL_DIRECTORIES} MATCHES ON )

# Install core files for the KratosMultiphysics python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${KRATOS_SOURCE_DIR}/kratos/python_interface/__init__.py" "KratosMultiphysics/__init__.py" )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${KRATOS_SOURCE_DIR}/kratos/python_interface/kratos_globals.py" "KratosMultiphysics/kratos_globals.py" )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${KRATOS_SOURCE_DIR}/kratos/python_interface/python_registry.py" "KratosMultiphysics/python_registry.py" )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${KRATOS_SOURCE_DIR}/kratos/python_interface/python_registry_lists.py" "KratosMultiphysics/python_registry_lists.py" )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${KRATOS_SOURCE_DIR}/kratos/python_interface/python_registry_utilities.py" "KratosMultiphysics/python_registry_utilities.py" )

# Install the libraries in the libs folder
install(FILES ${Boost_LIBRARIES} DESTINATION libs)
install(FILES ${EXTRA_INSTALL_LIBS} DESTINATION libs)

# Perform post-install tasks (add_subdirectory commands are guaranteed to be executed in order)
add_subdirectory("scripts/post_install")
