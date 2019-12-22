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

if(${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS} MATCHES ON)
    add_definitions(-DSTRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS)
    message("**** WARNING: Manually disabled advanced Constitutive laws ****")
    list(REMOVE_ITEM KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_advanced_constitutive_laws_to_python.cpp)
else(${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS} MATCHES ON)
    file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE
        ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE}
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_advanced_constitutive/*.cpp
    )
endif(${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS} MATCHES ON)

list(REMOVE_ITEM KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp)
list(INSERT KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE 0 ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp)

# StructuralMechanics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    if(NOT(DEFINED ${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS} OR ${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS}))
        file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES
            ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES}
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_advanced_constitutive_tests/*.cpp
    )
    endif(NOT(DEFINED ${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS} OR ${STRUCTURAL_DISABLE_ADVANCED_CONSTITUTIVE_LAWS}))
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosStructuralMechanicsCore SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosStructuralMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

# Define library Kratos which defines the basic python interface
pybind11_add_module(KratosStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStructuralMechanicsApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosStructuralMechanicsCore)
    cotire(KratosStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py KratosMultiphysics/StructuralMechanicsApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
)

# Install targets
install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )