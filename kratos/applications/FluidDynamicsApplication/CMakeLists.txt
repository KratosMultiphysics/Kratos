set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FLUID_DYNAMICS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_dynamics_python_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dynamic_vms.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spalart_allmaras_element.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/periodic_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFluidDynamicsApplication SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_SOURCES})
target_link_libraries(KratosFluidDynamicsApplication KratosCore )
set_target_properties(KratosFluidDynamicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")
