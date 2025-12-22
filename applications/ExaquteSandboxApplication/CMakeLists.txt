set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExaquteSandboxApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/MeshingApplication)

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/MeshingApplication )

# generate variables with the sources
set( KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application_variables.cpp

    # processes
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/weighted_divergence_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_divergencefree_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/calculate_divergence_process.cpp

    # utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/drag_and_moment_utilities.cpp
)

## define python interface sources
set( KRATOS_EXAQUTE_SANDBOX_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/exaqute_sandbox_python_application.cpp
)

# define exaqute sandbox application test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_EXAQUTE_SANDBOX_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(TARGET KratosExaquteSandboxCore SOURCES "${KRATOS_EXAQUTE_SANDBOX_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosExaquteSandboxCore SHARED ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES} ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_CORE})
target_link_libraries(KratosExaquteSandboxCore PRIVATE KratosCore KratosMeshingCore)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosExaquteSandboxApplication MODULE THIN_LTO ${KRATOS_EXAQUTE_SANDBOX_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosExaquteSandboxApplication PRIVATE KratosExaquteSandboxCore )
set_target_properties(KratosExaquteSandboxApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosExaquteSandboxCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ExaquteSandboxApplication.py KratosMultiphysics/ExaquteSandboxApplication/__init__.py )

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
install(TARGETS KratosExaquteSandboxCore DESTINATION libs )
install(TARGETS KratosExaquteSandboxApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosExaquteSandboxCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosExaquteSandboxApplication" PARENT_SCOPE)