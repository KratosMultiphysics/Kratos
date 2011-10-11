set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosTrilinosApplication ****")
 
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ )

#dependecies on other application -- TRILINOS specific
include_directories( ${CMAKE_SOURCE_DIR}/applications/incompressible_fluid_application )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )
 

#find mkl solver package -- this is our own find package!!
find_package(TRILINOS REQUIRED)
include_directories(${TRILINOS_INCLUDE_DIR})
link_directories(${TRILINOS_LIBRARIES})

##add definitions needed

## generate variables with the sources
set( KRATOS_TRILINOS_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/trilinos_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/trilinos_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_space_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_convergence_criterias_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_schemes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_trilinos_communicator_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_first.cpp
)


# ###############################################################
add_library(KratosTrilinosApplication SHARED ${KRATOS_TRILINOS_APPLICATION_SOURCES})
target_link_libraries(KratosTrilinosApplication KratosCore ${TRILINOS_LIBRARIES} )
set_target_properties(KratosTrilinosApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosTrilinosApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosTrilinosApplication DESTINATION libs )

