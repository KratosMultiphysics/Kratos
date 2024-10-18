kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/TrilinosApplication)

include_directories(
  ${KRATOS_SOURCE_DIR}/kratos
  ${KRATOS_SOURCE_DIR}/applications/RANSApplication
  ${KRATOS_SOURCE_DIR}/applications/TrilinosApplication
)

file(
  GLOB_RECURSE
  KRATOS_RANS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosRANSTrilinosExtension MODULE THIN_LTO ${KRATOS_RANS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )

target_include_directories(KratosRANSTrilinosExtension SYSTEM PUBLIC ${TRILINOS_INCLUDE_DIR})
target_link_libraries( KratosRANSTrilinosExtension PRIVATE KratosRANSCore KratosTrilinosCore KratosMPICore )
set_target_properties( KratosRANSTrilinosExtension PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosRANSTrilinosExtension PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosRANSTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosRANSTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosRANSTrilinosExtension" PARENT_SCOPE)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosRANSTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/RANSApplication")

set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosRANSTrilinosExtension" PARENT_SCOPE)