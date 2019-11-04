include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/applications/RANSModellingApplication
  ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication
)

if (NOT TRILINOS_FOUND)
  find_package(TRILINOS REQUIRED)
endif()

include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR})
link_directories(${TRILINOS_LIBRARIES})

file(
  GLOB_RECURSE
  KRATOS_RANS_MODELLING_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosRANSModellingTrilinosExtension MODULE THIN_LTO ${KRATOS_RANS_MODELLING_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosRANSModellingTrilinosExtension PRIVATE KratosRANSModellingCore KratosTrilinosCore KratosMPICore)
set_target_properties( KratosRANSModellingTrilinosExtension PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosRANSModellingTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosRANSModellingTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
if(USE_COTIRE MATCHES ON)
  cotire(KratosRANSModellingTrilinosExtension)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosRANSModellingTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/RANSModellingApplication")
