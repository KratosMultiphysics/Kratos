set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosTopologyOptimizationApplication ****")

################### PYBIND11
include(pybind11Tools)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/TopologyOptimizationApplication )

## TopologyOptimization Core sources
file(GLOB_RECURSE KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## TopologyOptimization testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## TopologyOptimization python interface sources
file(GLOB_RECURSE KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosTopologyOptimizationCore SHARED ${KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_CORE} ${KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosTopologyOptimizationCore PUBLIC KratosCore)
set_target_properties(KratosTopologyOptimizationCore PROPERTIES COMPILE_DEFINITIONS "TOPOLOGY_OPTIMIZATION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosTopologyOptimizationApplication MODULE THIN_LTO ${KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosTopologyOptimizationApplication PRIVATE KratosTopologyOptimizationCore)
set_target_properties(KratosTopologyOptimizationApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosTopologyOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosTopologyOptimizationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/TopologyOptimizationApplication.py KratosMultiphysics/TopologyOptimizationApplication/__init__.py )

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
install(TARGETS KratosTopologyOptimizationCore DESTINATION libs )
install(TARGETS KratosTopologyOptimizationApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosTopologyOptimizationCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosTopologyOptimizationApplication" PARENT_SCOPE)
