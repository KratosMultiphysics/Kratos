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
#        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spalart_allmaras_element.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spalart_allmaras.cpp

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_wall_condition.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/periodic_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFluidDynamicsApplication SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_SOURCES})
target_link_libraries(KratosFluidDynamicsApplication KratosCore )
set_target_properties(KratosFluidDynamicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsApplication.py" DESTINATION "${CMAKE_SOURCE_DIR}/KratosMultiphysics" )
