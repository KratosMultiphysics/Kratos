set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsHydraulicsApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/kratos/applications/FluidDynamicsApplication )

## FluidDynamicsHydraulics Core sources
file( GLOB_RECURSE KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_CORE_SOURCES

    # application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_hydraulics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_hydraulics_application_variables.cpp

    # utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp

)

## FluidDynamicsHydraulics python interface sources
file( GLOB_RECURSE KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## FluidDynamicsHydraulics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FluidDynamicsBiomedicalApplication core library (C++ parts)
add_library(KratosFluidDynamicsHydraulicsCore SHARED ${KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosFluidDynamicsHydraulicsCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties(KratosFluidDynamicsHydraulicsCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_HYDRAULICS_APPLICATION=EXPORT,API")

## define library Kratos which defines the basic python interface
pybind11_add_module(KratosFluidDynamicsHydraulicsApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_HYDRAULICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosFluidDynamicsHydraulicsApplication PRIVATE KratosFluidDynamicsHydraulicsCore)
set_target_properties(KratosFluidDynamicsHydraulicsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidDynamicsHydraulicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFluidDynamicsHydraulicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosFluidDynamicsHydraulicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosFluidDynamicsHydraulicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsHydraulicsApplication.py KratosMultiphysics/FluidDynamicsHydraulicsApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/FluidDynamicsHydraulicsApplication/python_registry_lists.py )

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
install(TARGETS KratosFluidDynamicsHydraulicsCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsHydraulicsApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidDynamicsHydraulicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidDynamicsHydraulicsApplication" PARENT_SCOPE)