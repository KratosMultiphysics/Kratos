set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES

  # application class and variables
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application_variables.cpp

  # utilities (compiled first because they are used within the element cpps)
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp

  # elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp

  # conditions
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp

  # fluid constitutive laws
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp

  # processes
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
)

## define python interface sources
file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## define fluid dynamics test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FluidDynamicsApplication core library (C++ parts)
add_library( KratosFluidDynamicsCore SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES} )
target_link_libraries( KratosFluidDynamicsCore PUBLIC KratosCore)
set_target_properties( KratosFluidDynamicsCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_APPLICATION=EXPORT,API")

## FluidDynamicsApplication python module
pybind11_add_module( KratosFluidDynamicsApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsApplication PRIVATE KratosFluidDynamicsCore)
set_target_properties( KratosFluidDynamicsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidDynamicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFluidDynamicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsApplication.py KratosMultiphysics/FluidDynamicsApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/FluidDynamicsApplication/python_registry_lists.py )

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
install(TARGETS KratosFluidDynamicsCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )

## If trilinos is avaliable, switch on the Fluid trilinos extension
if((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))
  message(STATUS "Enabling trilinos extension for KratosFluidDynamicsApplication")
  add_subdirectory(trilinos_extension)
endif((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidDynamicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidDynamicsApplication" PARENT_SCOPE)
