kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/TrilinosApplication)

include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication
  ${CMAKE_SOURCE_DIR}/applications/PoromechanicsApplication
  ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication
)

if (NOT TRILINOS_FOUND)
  find_package(TRILINOS REQUIRED)
endif()

include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR})
link_directories(${TRILINOS_LIBRARIES})

file(
  GLOB_RECURSE
  KRATOS_POROMECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosPoromechanicsTrilinosExtension MODULE THIN_LTO ${KRATOS_POROMECHANICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )
if(USE_COTIRE MATCHES ON)
  set_target_properties(KratosPoromechanicsTrilinosExtension PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY_UNITY")
  cotire(KratosPoromechanicsTrilinosExtension)
endif(USE_COTIRE MATCHES ON)
target_link_libraries( KratosPoromechanicsTrilinosExtension PRIVATE KratosPoromechanicsCore KratosStructuralMechanicsCore KratosTrilinosCore KratosMPICore ${TRILINOS_LIBRARIES})
set_target_properties( KratosPoromechanicsTrilinosExtension PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosPoromechanicsTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPoromechanicsTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
## installing the resulting libraries
install(TARGETS KratosPoromechanicsTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/PoromechanicsApplication")
