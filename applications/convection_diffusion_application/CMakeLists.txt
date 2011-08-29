set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConvectionDiffusionApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_CONVECTION_DIFFUSION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/convection_diffusion_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_convection_diffusion_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conv_diff_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conv_diff_change_of_phase_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conv_diff_3d.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_face2D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_face3D.cpp 
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosConvectionDiffusionApplication SHARED ${KRATOS_CONVECTION_DIFFUSION_APPLICATION_SOURCES})
target_link_libraries(KratosConvectionDiffusionApplication KratosCore )
set_target_properties(KratosConvectionDiffusionApplication PROPERTIES PREFIX "")
install(TARGETS KratosConvectionDiffusionApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")