set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSolidMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ConstitutiveModelsApplication )


#if(NOT (${CONSTITUTIVE_MODELS_APPLICATION} MATCHES "ON"))
#    message(SEND_ERROR "The ConsitutiveModels Application must be compiled before this one. Set CONSTITUTIVE_MODELS_APPLICATION=ON in the configure file.")
#endif(NOT (${CONSTITUTIVE_MODELS_APPLICATION} MATCHES "ON"))

## generate variables with the sources
set( KRATOS_SOLID_MECHANICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_cross_sections_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp

	# solid elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/linear_solid_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/large_displacement_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/small_displacement_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/small_displacement_bbar_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/axisymmetric_small_displacement_element.cpp
 	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/total_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/updated_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/axisymmetric_updated_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/large_displacement_U_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/updated_lagrangian_U_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/axisymmetric_updated_lagrangian_U_P_element.cpp

	# beam elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/small_displacement_beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/small_displacement_beam_element_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_emc_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_semc_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/geometrically_exact_rod_element.cpp

	# shell elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_elements/shell_thick_element_3D4N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_elements/shell_thin_element_3D3N.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/boundary_condition.cpp

	# load conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/line_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/surface_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/axisymmetric_point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/axisymmetric_line_load_condition.cpp

	# moment conditions		
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/point_moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/line_moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/surface_moment_condition.cpp

	# elastic conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/point_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/line_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/surface_elastic_condition.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/axisymmetric_point_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/axisymmetric_line_elastic_condition.cpp

	# utilites
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_cross_section.cpp
	
	## Constitutive Models Application ##  DEPRECATED CONSTITUTIVE LAWS IMPLEMENTATION
	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_U_P_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_U_P_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_orthotropic_3D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_axisym_2D_law.cpp


	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_J2_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_J2_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_J2_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_U_P_J2_axisym_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plastic_plane_stress_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_simo_ju_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_simo_ju_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_simo_ju_plane_stress_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_modified_mises_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_modified_mises_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/isotropic_damage_modified_mises_plane_stress_2D_law.cpp

	## Constitutive Models Application ##  DEPRECATED CONSTITUTIVE LAWS IMPLEMENTATION

	
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosSolidMechanicsApplication SHARED ${KRATOS_SOLID_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosSolidMechanicsApplication KratosCore)
##target_link_libraries(KratosSolidMechanicsApplication KratosCore KratosConstitutiveModelsApplication)
set_target_properties(KratosSolidMechanicsApplication PROPERTIES PREFIX "")
set_target_properties(KratosSolidMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "SOLID_MECHANICS_APPLICATION=EXPORT,API")
install(TARGETS KratosSolidMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosSolidMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/constitutive_laws DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/SolidMechanicsApplication.py" DESTINATION KratosMultiphysics )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/beam_profiles.csv" DESTINATION KratosMultiphysics )
