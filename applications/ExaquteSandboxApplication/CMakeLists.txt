set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExaquteSandboxApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

# generate variables with the sources
set( KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application_variables.cpp

    # processes
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/weighted_divergence_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/metrics_divergencefree_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/calculate_divergence_process.cpp
)

## define python interface sources
set( KRATOS_EXAQUTE_SANDBOX_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/exaqute_sandbox_python_application.cpp
)

# define exaqute sandbox application test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_EXAQUTE_SANDBOX_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosExaquteSandboxCore SHARED ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES} ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_CORE} ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosExaquteSandboxCore PRIVATE KratosCore KratosMeshingCore)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosExaquteSandboxApplication MODULE THIN_LTO ${KRATOS_EXAQUTE_SANDBOX_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosExaquteSandboxApplication PRIVATE KratosExaquteSandboxCore )
set_target_properties(KratosExaquteSandboxApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosExaquteSandboxApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ExaquteSandboxApplication.py KratosMultiphysics/ExaquteSandboxApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosExaquteSandboxCore DESTINATION libs )
install(TARGETS KratosExaquteSandboxApplication DESTINATION libs )
