include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/applications/ParticleMechanicsApplication
  ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication
)

if (NOT TRILINOS_FOUND)
  find_package(TRILINOS REQUIRED)
endif()

include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR})
link_directories(${TRILINOS_LIBRARIES})

file(
  GLOB_RECURSE
  KRATOS_PARTICLE_MECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosParticleMechanicsTrilinosExtension MODULE THIN_LTO ${KRATOS_PARTICLE_MECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )
if(USE_COTIRE MATCHES ON)
set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY_UNITY")
cotire(KratosParticleMechanicsTrilinosExtension)
endif(USE_COTIRE MATCHES ON)
target_link_libraries( KratosParticleMechanicsTrilinosExtension PRIVATE KratosParticleMechanicsCore KratosTrilinosCore KratosMPICore ${TRILINOS_LIBRARIES})
set_target_properties( KratosParticleMechanicsTrilinosExtension PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosParticleMechanicsTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
## installing the resulting libraries
install(TARGETS KratosParticleMechanicsTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/ParticleMechanicsApplication")