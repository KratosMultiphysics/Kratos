include_directories(
  ${KRATOS_SOURCE_DIR}/kratos
  ${KRATOS_SOURCE_DIR}/applications/CoSimulationApplication
)

file(
  GLOB_RECURSE
  KRATOS_CO_SIMULATION_MPI_EXTENSION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp
)

## Python module
pybind11_add_module( KratosCoSimulationMPIExtension MODULE THIN_LTO ${KRATOS_CO_SIMULATION_MPI_EXTENSION_PYTHON_INTERFACE_SOURCES} )

target_link_libraries( KratosCoSimulationMPIExtension PRIVATE KratosCoSimulationCore KratosMPICore co_sim_io_mpi )
set_target_properties( KratosCoSimulationMPIExtension PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosCoSimulationMPIExtension PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosCoSimulationMPIExtension PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosCoSimulationMPIExtension PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
## installing the resulting libraries
install(TARGETS KratosCoSimulationMPIExtension DESTINATION libs )
install(TARGETS co_sim_io_mpi DESTINATION libs )

###############################################################################
## install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MPIExtension.py" DESTINATION "KratosMultiphysics/CoSimulationApplication")

set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosCoSimulationMPIExtension" PARENT_SCOPE)
