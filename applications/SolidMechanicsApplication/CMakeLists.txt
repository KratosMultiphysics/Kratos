set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSolidMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_SOLID_MECHANICS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_conditions_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/isotropic_shell_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_3D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/large_displacement_3D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian_3D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian_2D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spatial_lagrangian_3D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spatial_lagrangian_2D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_3D_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_2D_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_force_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_force_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_moment_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_force_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_force_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_force_3D_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosSolidMechanicsApplication SHARED ${KRATOS_SOLID_MECHANICS_APPLICATION_SOURCES})
##target_link_libraries(KratosSolidMechanicsApplication KratosCore)
target_link_libraries(KratosSolidMechanicsApplication KratosCore)
##target_link_libraries(KratosSolidMechanicsApplication KratosCore triangle tetgen msuite tetMeshOpt)
set_target_properties(KratosSolidMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosSolidMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/SolidMechanicsApplication.py" DESTINATION KratosMultiphysics )
