set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRANSApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

MACRO(LIST_CONTAINS var value)
  SET(${var})
  FOREACH (value2 ${ARGN})
    IF (${value} STREQUAL ${value2})
      SET(${var} TRUE)
    ENDIF (${value} STREQUAL ${value2})
  ENDFOREACH (value2)
ENDMACRO(LIST_CONTAINS)

LIST_CONTAINS(uses_chimera "ChimeraApplication" ${LIST_OF_APPLICATIONS_TO_BE_COMPILED})
IF (uses_chimera)
    MESSAGE("RANSApplication is compiled with Chimera support.")
    kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ChimeraApplication)
ENDIF (uses_chimera)

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication )
IF (uses_chimera)
    include_directories( ${KRATOS_SOURCE_DIR}/applications/ChimeraApplication )
ENDIF (uses_chimera)

## RANSApplication Core instances sources
set(KRATOS_RANS_APPLICATION_CORE_INSTANCES
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_application_variables.cpp
)

## RANSApplication Core Utilities sources
file(GLOB_RECURSE KRATOS_RANS_APPLICATION_CORE_UTILITIES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## RANSApplication Core sources
file(GLOB_RECURSE KRATOS_RANS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
)

file(GLOB_RECURSE KRATOS_RANS_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## RANSApplication testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
file(GLOB_RECURSE KRATOS_RANS_APPLICATION_TESTING_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosRANSCore SHARED ${KRATOS_RANS_APPLICATION_CORE_INSTANCES} ${KRATOS_RANS_APPLICATION_CORE_UTILITIES} ${KRATOS_RANS_APPLICATION_CORE} ${KRATOS_RANS_APPLICATION_TESTING_SOURCES})

IF (uses_chimera)
    target_link_libraries(KratosRANSCore PUBLIC KratosCore KratosFluidDynamicsCore KratosChimeraApplicationCore)
else()
    target_link_libraries(KratosRANSCore PUBLIC KratosCore KratosFluidDynamicsCore)
ENDIF (uses_chimera)

set_target_properties(KratosRANSCore PROPERTIES COMPILE_DEFINITIONS "RANS_APPLICATION=EXPORT,API")

pybind11_add_module(KratosRANSApplication MODULE THIN_LTO ${KRATOS_RANS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosRANSApplication PRIVATE KratosRANSCore)
set_target_properties(KratosRANSApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosRANSCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosRANSApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRANSApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRANSApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/RANSApplication.py KratosMultiphysics/RANSApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
      )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosRANSCore DESTINATION libs )
install(TARGETS KratosRANSApplication DESTINATION libs )

if(${USE_MPI})
  message("Enabling trilinos extension for KratosRANSApplication")
  add_subdirectory(trilinos_extension)
endif(${USE_MPI})