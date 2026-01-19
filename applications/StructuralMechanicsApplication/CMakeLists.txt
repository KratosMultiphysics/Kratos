set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## StructuralMechanics Core sources
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constraints/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## StructuralMechanics python interface sources
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

list(REMOVE_ITEM KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp)
list(INSERT KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE 0 ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp)

# StructuralMechanics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(TARGET KratosStructuralMechanicsCore SOURCES "${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES}")
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosStructuralMechanicsCore SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE})
target_link_libraries(KratosStructuralMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

## StructuralMechanics benchmark sources. Disabled by default
if(${KRATOS_BUILD_BENCHMARK} MATCHES ON)
    file(GLOB_RECURSE KRATOS_BENCHMARK_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/*.cpp
    )

    foreach(file ${KRATOS_BENCHMARK_SOURCES})
        get_filename_component(filename ${file} NAME_WE)
        add_executable(${filename} ${file})
        target_link_libraries(${filename} PUBLIC KratosStructuralMechanicsCore benchmark::benchmark)
        set_target_properties(${filename} PROPERTIES COMPILE_DEFINITIONS "KRATOS_BENCHMARK=IMPORT,API")
        install(TARGETS ${filename} DESTINATION benchmark)
    endforeach(file ${KRATOS_BENCHMARK_SOURCES})
endif(${KRATOS_BUILD_BENCHMARK} MATCHES ON)

# Define library Kratos which defines the basic python interface
pybind11_add_module(KratosStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStructuralMechanicsApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosStructuralMechanicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py KratosMultiphysics/StructuralMechanicsApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosStructuralMechanicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosStructuralMechanicsApplication" PARENT_SCOPE)