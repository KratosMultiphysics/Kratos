set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsBiomedicalApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/kratos/applications/FluidDynamicsApplication)

## define application core sources (pure C++ parts)
file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_CORE_SOURCES

  # application class and variables
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_biomedical_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_biomedical_application_variables.cpp

  # utilities
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp

)

## define python interface sources
file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_BIOMEDICAL_PYTHON_INTERFACE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## define fluid dynamics biomedical test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FluidDynamicsBiomedicalApplication core library (C++ parts)
add_library( KratosFluidDynamicsBiomedicalCore SHARED ${KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES} )
target_link_libraries( KratosFluidDynamicsBiomedicalCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties( KratosFluidDynamicsBiomedicalCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_BIOMEDICAL_APPLICATION=EXPORT,API")

## FluidDynamicsBiomedicalApplication python module
pybind11_add_module( KratosFluidDynamicsBiomedicalApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_BIOMEDICAL_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsBiomedicalApplication PRIVATE KratosFluidDynamicsBiomedicalCore)
set_target_properties( KratosFluidDynamicsBiomedicalApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidDynamicsBiomedicalCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFluidDynamicsBiomedicalApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsBiomedicalApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsBiomedicalApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsBiomedicalApplication.py KratosMultiphysics/FluidDynamicsBiomedicalApplication/__init__.py )

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
install(TARGETS KratosFluidDynamicsBiomedicalCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsBiomedicalApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidDynamicsBiomedicalCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidDynamicsBiomedicalApplication" PARENT_SCOPE)
