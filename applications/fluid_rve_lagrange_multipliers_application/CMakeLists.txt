set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidRveLagrangeMultipliersApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FLUIDRVELAGRANGEMULTIPLIERS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//fluid_rve_lagrange_multipliers_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_rve_lagrange_multipliers_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_mean_velocity_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_velocity_gradients_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/inverse_tangent_velocity_periodic_condition_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/inverse_normal_velocity_periodic_condition_2d.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFluidRveLagrangeMultipliersApplication SHARED ${KRATOS_FLUIDRVELAGRANGEMULTIPLIERS_APPLICATION_SOURCES})
target_link_libraries(KratosFluidRveLagrangeMultipliersApplication KratosCore )
set_target_properties(KratosFluidRveLagrangeMultipliersApplication PROPERTIES PREFIX "")
install(TARGETS KratosFluidRveLagrangeMultipliersApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFluidRveLagrangeMultipliersApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidRveLagrangeMultipliersApplication.py" DESTINATION KratosMultiphysics )
