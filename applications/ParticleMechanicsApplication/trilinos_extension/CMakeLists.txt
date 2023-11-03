kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/TrilinosApplication)

include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/applications/ParticleMechanicsApplication
  ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication
)

file(
  GLOB_RECURSE
  KRATOS_PARTICLE_MECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_builder_and_solver/*.cpp
)

## Python module
pybind11_add_module( KratosParticleMechanicsTrilinosExtension MODULE THIN_LTO ${KRATOS_PARTICLE_MECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )

target_link_libraries( KratosParticleMechanicsTrilinosExtension PRIVATE KratosParticleMechanicsCore KratosTrilinosCore KratosMPICore )
set_target_properties( KratosParticleMechanicsTrilinosExtension PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosParticleMechanicsTrilinosExtension" PARENT_SCOPE)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosParticleMechanicsTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/ParticleMechanicsApplication")

set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosParticleMechanicsTrilinosExtension" PARENT_SCOPE)
