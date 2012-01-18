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

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosConvectionDiffusionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosConvectionDiffusionApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ConvectionDiffusionApplication.py" DESTINATION KratosMultiphysics )
