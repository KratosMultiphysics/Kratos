set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosKaHIPApplication ****")

include(pybind11Tools)

# ── Locate KaHIP ──────────────────────────────────────────────────────────────
# FindKaHIP.cmake lives in cmake_modules/ (already on CMAKE_MODULE_PATH).
# Search order: KAHIP_ROOT cache/env → external_libraries/KaHIP/deploy (auto-built)
# → external_libraries/KaHIP/build → system paths → pkg-config.

option(KAHIP_DOWNLOAD_IF_NOT_FOUND
    "Automatically clone and build KaHIP into external_libraries/KaHIP if not found" ON)

find_package(KaHIP QUIET)

if(NOT KAHIP_FOUND)
    if(NOT KAHIP_DOWNLOAD_IF_NOT_FOUND)
        message(FATAL_ERROR
            "KaHIP not found and KAHIP_DOWNLOAD_IF_NOT_FOUND is OFF.\n"
            "Set KAHIP_ROOT to a KaHIP build/install directory, or enable "
            "KAHIP_DOWNLOAD_IF_NOT_FOUND to auto-build it.")
    endif()

    # ── Auto-download and build KaHIP ─────────────────────────────────────────
    set(_kahip_src_dir   "${KRATOS_SOURCE_DIR}/external_libraries/KaHIP")
    set(_kahip_build_dir "${_kahip_src_dir}/build")
    set(_kahip_deploy_dir "${_kahip_src_dir}/deploy")

    # Clone KaHIP source if not already present
    if(NOT EXISTS "${_kahip_src_dir}/CMakeLists.txt")
        message(STATUS "KaHIP not found – cloning from GitHub into ${_kahip_src_dir} ...")
        execute_process(
            COMMAND git clone --depth 1
                    https://github.com/KaHIP/KaHIP.git "${_kahip_src_dir}"
            RESULT_VARIABLE _kahip_clone_rc
            ERROR_VARIABLE  _kahip_clone_err)
        if(NOT _kahip_clone_rc EQUAL 0)
            message(FATAL_ERROR "Failed to clone KaHIP:\n${_kahip_clone_err}")
        endif()
        message(STATUS "KaHIP cloned successfully.")
    else()
        message(STATUS "KaHIP source already present at ${_kahip_src_dir}.")
    endif()

    # Copy the provided configure.sh into the KaHIP build directory so that
    # SOURCE_DIR (parent of SCRIPT_DIR) resolves to the KaHIP source root.
    file(MAKE_DIRECTORY "${_kahip_build_dir}")
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/KaHIP/configure.sh"
        "${_kahip_build_dir}/configure.sh"
        COPYONLY)

    # Forward MPI availability to the KaHIP build
    if(USE_MPI MATCHES ON AND MPI_FOUND)
        set(_kahip_mpi_args "")
    else()
        set(_kahip_mpi_args "--no-mpi")
    endif()

    # Build KaHIP at CMake configure time (synchronous)
    message(STATUS "Building KaHIP – this may take a few minutes ...")
    execute_process(
        COMMAND bash configure.sh --no-native ${_kahip_mpi_args}
        WORKING_DIRECTORY "${_kahip_build_dir}"
        RESULT_VARIABLE _kahip_build_rc
        OUTPUT_VARIABLE _kahip_build_out
        ERROR_VARIABLE  _kahip_build_err)
    if(NOT _kahip_build_rc EQUAL 0)
        message(FATAL_ERROR
            "KaHIP build failed (exit code ${_kahip_build_rc}):\n${_kahip_build_err}")
    endif()
    message(STATUS "KaHIP built successfully. Artifacts in ${_kahip_deploy_dir}")

    # configure.sh deploys to <source>/deploy by default; point KAHIP_ROOT there
    set(KAHIP_ROOT "${_kahip_deploy_dir}" CACHE PATH "KaHIP deploy directory" FORCE)

    # Re-run find now that the library is in place
    find_package(KaHIP REQUIRED)
endif()

message(STATUS "KaHIP include dir : ${KAHIP_INCLUDE_DIRS}")
message(STATUS "KaHIP libraries   : ${KAHIP_LIBRARIES}")

# ── Optional ParHIP (MPI-parallel partitioning) ───────────────────────────────
option(KAHIP_USE_PARHIP "Enable ParHIP for MPI-parallel mesh partitioning" ON)

if(KAHIP_USE_PARHIP AND KAHIP_PARHIP_LIBRARIES)
    message(STATUS "ParHIP libraries  : ${KAHIP_PARHIP_LIBRARIES}")
    set(_KAHIP_HAS_PARHIP TRUE)
else()
    if(KAHIP_USE_PARHIP AND NOT KAHIP_PARHIP_LIBRARIES)
        message(WARNING "KAHIP_USE_PARHIP is ON but parhip_interface library was not found. "
                        "MPI-parallel partitioning will not be available.")
    endif()
    set(_KAHIP_HAS_PARHIP FALSE)
endif()

# ── Optional 64-bit index mode ────────────────────────────────────────────────
option(KAHIP_64BIT "Build KaHIPApplication with 64-bit kahip_idx (requires KaHIP built with -D64BITMODE=ON)" OFF)

# ── Include directories ───────────────────────────────────────────────────────
include_directories(${KRATOS_SOURCE_DIR}/kratos)
include_directories(${KAHIP_INCLUDE_DIRS})

# ── Core source files ─────────────────────────────────────────────────────────
file(GLOB KRATOS_KAHIP_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/kahip_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/kahip_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_modeler/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

# ── Python interface sources ──────────────────────────────────────────────────
file(GLOB_RECURSE KRATOS_KAHIP_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# ── C++ unit tests ────────────────────────────────────────────────────────────
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB KRATOS_KAHIP_APPLICATION_TESTING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(
        TARGET  KratosKaHIPCore
        SOURCES "${KRATOS_KAHIP_APPLICATION_TESTING_SOURCES}"
        USE_CUSTOM_MAIN)
endif()

# ── Core shared library ───────────────────────────────────────────────────────
add_library(KratosKaHIPCore SHARED ${KRATOS_KAHIP_APPLICATION_CORE})
target_link_libraries(KratosKaHIPCore PUBLIC KratosCore ${KAHIP_LIBRARIES})

# Link ParHIP if available
if(_KAHIP_HAS_PARHIP)
    target_link_libraries(KratosKaHIPCore PUBLIC ${KAHIP_PARHIP_LIBRARIES})
    target_compile_definitions(KratosKaHIPCore PUBLIC KAHIP_USE_PARHIP)
endif()

# 64-bit index mode
if(KAHIP_64BIT)
    target_compile_definitions(KratosKaHIPCore PUBLIC KAHIP_64BIT)
endif()

# Compiler export macro
set_target_properties(KratosKaHIPCore PROPERTIES
    COMPILE_DEFINITIONS "KAHIP_APPLICATION=EXPORT,API")

# ── Python extension module ───────────────────────────────────────────────────
pybind11_add_module(KratosKaHIPApplication MODULE THIN_LTO
    ${KRATOS_KAHIP_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosKaHIPApplication PRIVATE KratosKaHIPCore)
set_target_properties(KratosKaHIPApplication PROPERTIES PREFIX "")

# Unity build batch size
if(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosKaHIPCore PROPERTIES
        UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosKaHIPApplication PROPERTIES
        UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
endif()

# Platform-specific suffix adjustments
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosKaHIPApplication PROPERTIES SUFFIX .pyd)
endif()
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosKaHIPApplication PROPERTIES SUFFIX .so)
endif()

# ── Python package install ────────────────────────────────────────────────────
kratos_python_install(
    ${INSTALL_PYTHON_USING_LINKS}
    ${CMAKE_CURRENT_SOURCE_DIR}/KaHIPApplication.py
    KratosMultiphysics/KaHIPApplication/__init__.py)

kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS}
    "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/"
    KratosMultiphysics/KaHIPApplication)

# ── Test data install ─────────────────────────────────────────────────────────
if(${INSTALL_TESTING_FILES} MATCHES ON)
    get_filename_component(CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests
        DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git"  EXCLUDE
        PATTERN "*.c"    EXCLUDE
        PATTERN "*.h"    EXCLUDE
        PATTERN "*.cpp"  EXCLUDE
        PATTERN "*.hpp"  EXCLUDE)
endif()

# ── Library install ───────────────────────────────────────────────────────────
install(TARGETS KratosKaHIPCore DESTINATION libs)
install(TARGETS KratosKaHIPApplication DESTINATION libs)

# Install KaHIP shared libraries into the Kratos libs folder so the runtime
# linker finds them without any LD_LIBRARY_PATH extension beyond the Kratos one.
get_filename_component(KAHIP_LIB_DIR "${KAHIP_LIBRARIES}" DIRECTORY)
file(GLOB KAHIP_SHARED_LIBS "${KAHIP_LIB_DIR}/libkahip*.so*")
foreach(LIB_FILE ${KAHIP_SHARED_LIBS})
    install(FILES "${LIB_FILE}" DESTINATION libs)
endforeach()

if(_KAHIP_HAS_PARHIP)
    install(FILES "${KAHIP_PARHIP_LIBRARIES}" DESTINATION libs)
endif()

# ── Register with Kratos build system ────────────────────────────────────────
set(KRATOS_KERNEL
    "${KRATOS_KERNEL};KratosKaHIPCore"
    PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE
    "${KRATOS_PYTHON_INTERFACE};KratosKaHIPApplication"
    PARENT_SCOPE)
