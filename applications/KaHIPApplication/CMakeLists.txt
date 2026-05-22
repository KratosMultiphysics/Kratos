set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosKaHIPApplication ****")

include(pybind11Tools)

# ── Locate KaHIP ──────────────────────────────────────────────────────────────
# Extend the CMake module search path to find our bundled FindKaHIP.cmake first.
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Default: prefer the in-tree KaHIP build shipped under external_libraries/KaHIP
if(NOT DEFINED KAHIP_ROOT AND NOT DEFINED ENV{KAHIP_ROOT})
    set(KAHIP_ROOT "${KRATOS_SOURCE_DIR}/external_libraries/KaHIP/build"
        CACHE PATH "KaHIP installation or build directory")
endif()

find_package(KaHIP REQUIRED)

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

# Install the KaHIP shared library so it can be found at runtime
get_filename_component(KAHIP_LIB_DIR "${KAHIP_LIBRARIES}" DIRECTORY)
file(GLOB KAHIP_SHARED_LIBS "${KAHIP_LIB_DIR}/libkahip*.so*")
foreach(LIB_FILE ${KAHIP_SHARED_LIBS})
    install(FILES "${LIB_FILE}" DESTINATION libs)
endforeach()

# ── Register with Kratos build system ────────────────────────────────────────
set(KRATOS_KERNEL
    "${KRATOS_KERNEL};KratosKaHIPCore"
    PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE
    "${KRATOS_PYTHON_INTERFACE};KratosKaHIPApplication"
    PARENT_SCOPE)
