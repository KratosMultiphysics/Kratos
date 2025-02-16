cmake_minimum_required(VERSION 3.14)
# CMake version compatibility
# TODO: Remove when bumping cmake >= 3.25
if (POLICY CMP0140)
    # Enables: return(PROPAGATE)
    cmake_policy(SET CMP0140 NEW)
endif ()

#[==============================================================================================[
#                                    Basic project definition                                   #
]==============================================================================================]

# TODO: Version 3, rename the project and namespace to something more compact
project(nlohmann_json_schema_validator
        VERSION 2.3.0
        DESCRIPTION "Json validator for nlohmann::json library"
        HOMEPAGE_URL "https://github.com/pboettch/json-schema-validator"
        LANGUAGES CXX)
# TODO: Remove when bumping cmake >= 3.21
if (NOT DEFINED nlohmann_json_schema_validator_IS_TOP_LEVEL)
    if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        set(PROJECT_IS_TOP_LEVEL ON)
    else ()
        set(PROJECT_IS_TOP_LEVEL OFF)
    endif ()
endif ()

#[==============================================================================================[
#                                            Options                                            #
]==============================================================================================]

option(JSON_VALIDATOR_INSTALL "JsonValidator: Install targets" ${PROJECT_IS_TOP_LEVEL})
option(JSON_VALIDATOR_BUILD_TESTS "JsonValidator: Build tests" ${PROJECT_IS_TOP_LEVEL})
option(JSON_VALIDATOR_BUILD_EXAMPLES "JsonValidator: Build examples" ${PROJECT_IS_TOP_LEVEL})
option(JSON_VALIDATOR_SHARED_LIBS "JsonValidator: Build as shared library" ${PROJECT_IS_TOP_LEVEL})
option(JSON_VALIDATOR_TEST_COVERAGE "JsonValidator: Build with test coverage" OFF)
mark_as_advanced(JSON_VALIDATOR_TEST_COVERAGE)
# Get a default JSON_FETCH_VERSION from environment variables to workaround the CI
if (DEFINED ENV{NLOHMANN_JSON_VERSION})
    set(JSON_FETCH_VERSION_DEFAULT $ENV{NLOHMANN_JSON_VERSION})
else ()
    set(JSON_FETCH_VERSION_DEFAULT v3.11.2)
endif ()
set(JSON_FETCH_VERSION ${JSON_FETCH_VERSION_DEFAULT} CACHE STRING "Fetch nlohmann::json version")

#[==============================================================================================[
#                                     Project configuration                                     #
]==============================================================================================]

# Include cmake modules
include(FetchContent)
if (JSON_VALIDATOR_INSTALL)
    include(GNUInstallDirs)
    include(CMakePackageConfigHelpers)
endif ()

# Default to release build
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

# Enable cmake's BUILD_SHARED_LIBS
set(BUILD_SHARED_LIBS ${nlohmann_json_schema_validator_SHARED_LIBS})

if (JSON_VALIDATOR_TEST_COVERAGE)
    if (CMAKE_CXX_COMPILER_ID STREQUAL Clang)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL GNU)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    else ()
        message(WARNING
                "JsonValidator: Other toolchain coverage flags unknown.\n"
                "Using --coverage as default")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
    endif ()
endif ()

#[==============================================================================================[
#                                       External packages                                       #
]==============================================================================================]

set(fetch_packages "")
if (NOT TARGET nlohmann_json)
  # Fetch/Find nlohmann_json
  # TODO: Remove when bumping cmake >= 3.24
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.24)
      FetchContent_Declare(nlohmann_json
              GIT_REPOSITORY https://github.com/nlohmann/json
              GIT_TAG ${JSON_FETCH_VERSION}
              FIND_PACKAGE_ARGS
              )
      list(APPEND fetch_packages nlohmann_json)
  else ()
      # Try to get system installed version
      find_package(nlohmann_json QUIET)
      if (NOT nlohmann_json_FOUND)
          # If failed fetch the desired version
          FetchContent_Declare(nlohmann_json
                  GIT_REPOSITORY https://github.com/nlohmann/json
                  GIT_TAG ${JSON_FETCH_VERSION}
                  )
          list(APPEND fetch_packages nlohmann_json)
      endif ()
  endif ()
endif ()

# Handle configure flags
if (JSON_VALIDATOR_INSTALL)
    # TODO: This is not ideal, this package should not be installing nlohmann::json
    # Currently required in order to satisfy cmake exporter
    set(JSON_Install ON CACHE BOOL "")
endif ()

# Get all dependencies
FetchContent_MakeAvailable(${fetch_packages})
if (JSON_VALIDATOR_INSTALL AND NOT nlohmann_json_FOUND AND JSON_Install)
    # TODO: This is not ideal
    message(WARNING
            "JsonValidator: No nlohmann::json found on the system and nlohmann_json_schema_validator will be installed\n"
            "This will also install nlohmann::json in its typical installation path\n"
            "This is not ideal because it might overwrite system installed")
endif ()

#[==============================================================================================[
#                                        Main definition                                        #
]==============================================================================================]

message(STATUS "JsonValidator: Configured for ${CMAKE_BUILD_TYPE}")
if (DEFINED nlohmann_json_VERSION)
    message(STATUS "JsonValidator: Using nlohmann/json version: ${nlohmann_json_VERSION}")
else ()
    message(STATUS "JsonValidator: nlohmann_json_VERSION is not set. Possible value: ${JSON_FETCH_VERSION}")
endif ()

## Main targets
add_library(nlohmann_json_schema_validator)
add_library(nlohmann_json_schema_validator::validator ALIAS nlohmann_json_schema_validator)
set_target_properties(nlohmann_json_schema_validator PROPERTIES
        VERSION ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
        EXPORT_NAME validator
        # TODO: Version 3, simplify the library name
#        OUTPUT_NAME nlohmann_json_validator
        )

# Main definitions in here
add_subdirectory(src)

# Enable examples

# Enable testings
if (JSON_VALIDATOR_BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif ()

if (JSON_VALIDATOR_BUILD_EXAMPLES)
    add_subdirectory(example)
endif ()


#[==============================================================================================[
#                                       Install or Export                                       #
]==============================================================================================]

if (JSON_VALIDATOR_INSTALL)
    # Note other install targets found in subdirectories
    # Here mostly the cmake boilerplate are set
    write_basic_package_version_file(nlohmann_json_schema_validatorConfigVersion.cmake
            VERSION ${PROJECT_VERSION}
            COMPATIBILITY SameMajorVersion
            )
    configure_package_config_file(cmake/nlohmann_json_schema_validatorConfig.cmake.in
            nlohmann_json_schema_validatorConfig.cmake
            INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nlohmann_json_schema_validator
            )

    # Install Targets files
    export(EXPORT nlohmann_json_schema_validatorTargets
            NAMESPACE nlohmann_json_schema_validator::
            FILE nlohmann_json_schema_validatorTargets.cmake
            )
    install(EXPORT nlohmann_json_schema_validatorTargets
            FILE nlohmann_json_schema_validatorTargets.cmake
            NAMESPACE nlohmann_json_schema_validator::
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nlohmann_json_schema_validator
            COMPONENT nlohmann_json_schema_validator_Development
            )
    # Install cmake export files
    install(FILES
            ${CMAKE_CURRENT_BINARY_DIR}/nlohmann_json_schema_validatorConfig.cmake
            ${CMAKE_CURRENT_BINARY_DIR}/nlohmann_json_schema_validatorConfigVersion.cmake
            DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/nlohmann_json_schema_validator
            COMPONENT nlohmann_json_schema_validator_Development
            )
endif ()

# Handle the project being included externally (e.g. FetchContent)
if (NOT PROJECT_IS_TOP_LEVEL)
    # Export variables set in nlohmann_json_schema_validatorConfig.cmake
    # TODO: Remove when bumping cmake >= 3.25
    if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.25)
        return(PROPAGATE
                nlohmann_json_schema_validator_VERSION
                nlohmann_json_schema_validator_VERSION_MAJOR
                nlohmann_json_schema_validator_VERSION_MINOR
                nlohmann_json_schema_validator_VERSION_PATCH
                nlohmann_json_schema_validator_VERSION_TWEAK
                )
    else ()
        set(nlohmann_json_schema_validator_VERSION ${nlohmann_json_schema_validator_VERSION} PARENT_SCOPE)
        set(nlohmann_json_schema_validator_VERSION_MAJOR ${nlohmann_json_schema_validator_VERSION_MAJOR} PARENT_SCOPE)
        set(nlohmann_json_schema_validator_VERSION_MINOR ${nlohmann_json_schema_validator_VERSION_MINOR} PARENT_SCOPE)
        set(nlohmann_json_schema_validator_VERSION_PATCH ${nlohmann_json_schema_validator_VERSION_PATCH} PARENT_SCOPE)
        set(nlohmann_json_schema_validator_VERSION_TWEAK ${nlohmann_json_schema_validator_VERSION_TWEAK} PARENT_SCOPE)
    endif ()
endif ()
