set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication ) 
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication ) 

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//structural_mechanics_application.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}//structural_mechanics_application_variables.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/structural_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_cross_sections_to_python.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_cross_section.cpp
        
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_beam_element_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/membrane_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/isotropic_shell_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_element_3D4N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thin_element_3D3N.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_moment_3D_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosStructuralMechanicsApplication SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosStructuralMechanicsApplication KratosCore KratosSolidMechanicsApplication)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/beam_profiles.csv" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/beam_profiles.csv" DESTINATION applications/${CURRENT_DIR_NAME}/python_scripts )
