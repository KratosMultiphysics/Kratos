set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSolidMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ConstitutiveModelsApplication)

## Generate variables with the sources
set( KRATOS_SOLID_MECHANICS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/solid_mechanics_application_variables.cpp

	# Time integration methods
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/solution_local_flags.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/static_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/newmark_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/simo_step_rotation_method.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_solvers/time_integration_methods/emc_step_rotation_method.cpp

	# Solid elements
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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/large_displacement_V_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/updated_lagrangian_V_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/large_displacement_segregated_V_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_elements/updated_lagrangian_segregated_V_P_element.cpp

	# Beam elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/small_displacement_beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/small_displacement_beam_element_3D2N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_emc_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/large_displacement_beam_semc_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_elements/geometrically_exact_rod_element.cpp

	# Shell elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_elements/shell_thick_element_3D4N.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_elements/shell_thin_element_3D3N.cpp

	# Thermal elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_elements/thermal_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_elements/axisymmetric_thermal_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/boundary_condition.cpp

	# Load conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/line_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/surface_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/axisymmetric_point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_conditions/axisymmetric_line_load_condition.cpp

	# Moment conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/point_moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/line_moment_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/moment_conditions/surface_moment_condition.cpp

	# Elastic conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/point_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/line_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/surface_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/axisymmetric_point_elastic_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/elastic_conditions/axisymmetric_line_elastic_condition.cpp

	# Thermal conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_conditions/line_heat_flux_condition.cpp

	# Utilites
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_cross_section.cpp

	# Processes
	${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/assign_scalar_variable_to_entities_process.cpp


	##  DEPRECATED CONSTITUTIVE LAWS IMPLEMENTATION ##
	# Constitutive laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_thermal_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/johnson_cook_thermal_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/baker_johnson_cook_thermal_hardening_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/mises_huber_thermal_yield_criterion.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_linear_rate_dependent_plastic_flow_rule.cpp

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

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_johnson_cook_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_baker_johnson_cook_plane_strain_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_johnson_cook_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_johnson_cook_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_thermal_U_P_baker_johnson_cook_plane_strain_2D_law.cpp

	## DEPRECATED CONSTITUTIVE LAWS IMPLEMENTATION ##
)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_SOLID_MECHANICS_APPLICATION_TESTS ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)


## Generate variables with the sources
set( KRATOS_SOLID_MECHANICS_APPLICATION_PYTHON_INTERFACE
        # Custom python
  	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_cross_sections_to_python.cpp
)

add_library(KratosSolidMechanicsCore SHARED ${KRATOS_SOLID_MECHANICS_APPLICATION_CORE} ${KRATOS_SOLID_MECHANICS_APPLICATION_TESTS})
#target_link_libraries(KratosSolidMechanicsCore PUBLIC KratosCore KratosConstitutiveModelsCore)
target_link_libraries(KratosSolidMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosSolidMechanicsCore PROPERTIES COMPILE_DEFINITIONS "SOLID_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosSolidMechanicsApplication MODULE THIN_LTO ${KRATOS_SOLID_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosSolidMechanicsApplication PRIVATE KratosSolidMechanicsCore)
set_target_properties(KratosSolidMechanicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosSolidMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosSolidMechanicsCore)
    cotire(KratosSolidMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/SolidMechanicsApplication.py KratosMultiphysics/SolidMechanicsApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosSolidMechanicsCore DESTINATION libs )
install(TARGETS KratosSolidMechanicsApplication DESTINATION libs )