include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/applications/ChimeraApplication
  ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication
)

if (NOT TRILINOS_FOUND)
  find_package(TRILINOS REQUIRED)
endif()

include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR})
link_directories(${TRILINOS_LIBRARIES})

file(
  GLOB_RECURSE
  KRATOS_CHIMERA_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

file(
  GLOB_RECURSE
  TRILINOS_TEST_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/trilinos_test.cpp
)


## Python module
pybind11_add_module( KratosChimeraTrilinosExtension MODULE THIN_LTO ${KRATOS_CHIMERA_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )
if(USE_COTIRE MATCHES ON)
set_target_properties(KratosChimeraTrilinosExtension PROPERTIES COTIRE_UNITY_LINK_LIBRARIES_INIT "COPY_UNITY")
cotire(KratosChimeraTrilinosExtension)
endif(USE_COTIRE MATCHES ON)
target_link_libraries( KratosChimeraTrilinosExtension PRIVATE KratosChimeraApplicationCore KratosTrilinosCore KratosMPICore ${TRILINOS_LIBRARIES})
set_target_properties( KratosChimeraTrilinosExtension PROPERTIES PREFIX "")

add_executable (trilinos_test ${TRILINOS_TEST_SOURCES})
target_link_libraries (trilinos_test PRIVATE KratosTrilinosCore KratosMPICore ${TRILINOS_LIBRARIES})
###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosChimeraTrilinosExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosChimeraTrilinosExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
## installing the resulting libraries
install(TARGETS KratosChimeraTrilinosExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TrilinosExtension.py" DESTINATION "KratosMultiphysics/ChimeraApplication")
