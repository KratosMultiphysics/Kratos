set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosThermoMechanicalApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_THERMO_MECHANICAL_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/thermo_mechanical_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/thermo_mechanical_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp

# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/proj_dirichlet_cond.cpp		
# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fluid3d_neumann.cpp
# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/no_slip_condition_2d.cpp
# 	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic2d_neumann.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/heat_contact_2d.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/heat_contact_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/SUPG_conv_diff_2d.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/SUPG_conv_diff_3d.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/SUPG_conv_3d.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_face2d.cpp  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_face3d.cpp  

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosThermoMechanicalApplication SHARED ${KRATOS_THERMO_MECHANICAL_APPLICATION_SOURCES})
target_link_libraries(KratosThermoMechanicalApplication KratosCore )
set_target_properties(KratosThermoMechanicalApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosThermoMechanicalApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosThermoMechanicalApplication DESTINATION libs )


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)



# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ThermoMechanicalApplication.py" DESTINATION KratosMultiphysics )
