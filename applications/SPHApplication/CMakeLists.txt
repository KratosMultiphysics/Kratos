set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSPHApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## SPH Core sources
file(GLOB_RECURSE KRATOS_SPH_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/sph_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sph_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/custom_kernels/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
)

## SPH testing sources
#if(${KRATOS_BUILD_TESTING} MATCHES ON)
#    file(GLOB_RECURSE KRATOS_SPH_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
#endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## SPH python interface sources
file(GLOB_RECURSE KRATOS_SPH_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosSPHCore SHARED ${KRATOS_SPH_APPLICATION_CORE} ${KRATOS_SPH_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosSPHCore PUBLIC KratosCore)
set_target_properties(KratosSPHCore PROPERTIES COMPILE_DEFINITIONS "SPH_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosSPHApplication MODULE THIN_LTO ${KRATOS_SPH_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosSPHApplication PRIVATE KratosSPHCore)
set_target_properties(KratosSPHApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosSPHApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosSPHApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/SPHApplication.py KratosMultiphysics/SPHApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/SPHApplication/python_registry_lists.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
# Test files removed 

# Install targets
install(TARGETS KratosSPHCore DESTINATION libs )
install(TARGETS KratosSPHApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosSPHCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosSPHApplication" PARENT_SCOPE)