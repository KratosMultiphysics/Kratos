kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/TrilinosApplication)

include_directories(
  ${KRATOS_SOURCE_DIR}/kratos
  ${KRATOS_SOURCE_DIR}/applications/MappingApplication
  ${KRATOS_SOURCE_DIR}/applications/TrilinosApplication
)

file(
  GLOB_RECURSE
  KRATOS_MAPPING_MPI_EXTENSION_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_searching/*.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## define test sources (cpp tests)
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB KRATOS_MAPPING_MPI_EXTENSION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cpp_tests/*.cpp)
  kratos_add_gtests(TARGET KratosMappingMPICore SOURCES "${KRATOS_MAPPING_MPI_EXTENSION_TESTING_SOURCES}" USE_MPI)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

file(
  GLOB_RECURSE
  KRATOS_MAPPING_MPI_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

add_library(KratosMappingMPICore SHARED ${KRATOS_MAPPING_MPI_EXTENSION_CORE_SOURCES})
target_include_directories(KratosMappingMPICore SYSTEM PUBLIC ${TRILINOS_INCLUDE_DIR})
target_link_libraries(KratosMappingMPICore PUBLIC KratosMPICore KratosCore KratosMappingCore KratosTrilinosCore ${MPI_LIBRARIES} ${TRILINOS_LIBRARIES} )
set_target_properties(KratosMappingMPICore PROPERTIES COMPILE_DEFINITIONS "MAPPING_MPI_APPLICATION=EXPORT,API")

## Python module
pybind11_add_module( KratosMappingMPIExtension MODULE THIN_LTO ${KRATOS_MAPPING_MPI_EXTENSION_PYTHON_INTERFACE_SOURCES} )

target_link_libraries( KratosMappingMPIExtension PRIVATE KratosMappingMPICore )
set_target_properties( KratosMappingMPIExtension PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMappingMPIExtension PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosMappingMPIExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMappingMPIExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMappingMPICore" PARENT_SCOPE)
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMappingMPIExtension" PARENT_SCOPE)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosMappingMPICore DESTINATION libs )
install(TARGETS KratosMappingMPIExtension DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MPIExtension.py" DESTINATION "KratosMultiphysics/MappingApplication")

set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMappingMPIExtension" PARENT_SCOPE)
