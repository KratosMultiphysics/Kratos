kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/TrilinosApplication)

include_directories(
  ${KRATOS_SOURCE_DIR}/kratos
  ${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication
  ${KRATOS_SOURCE_DIR}/applications/TrilinosApplication
)

file(
  GLOB_RECURSE
  KRATOS_FLUID_DYNAMICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosFluidDynamicsTrilinosExtension MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )

target_include_directories(KratosFluidDynamicsTrilinosExtension SYSTEM PUBLIC ${TRILINOS_INCLUDE_DIR})
target_link_libraries( KratosFluidDynamicsTrilinosExtension PRIVATE KratosFluidDynamicsCore KratosTrilinosCore KratosMPICore )
set_target_properties( KratosFluidDynamicsTrilinosExtension PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidDynamicsTrilinosExtension PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidDynamicsTrilinosExtension" PARENT_SCOPE)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosFluidDynamicsTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/FluidDynamicsApplication")

set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidDynamicsTrilinosExtension" PARENT_SCOPE)