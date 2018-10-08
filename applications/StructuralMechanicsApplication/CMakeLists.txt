set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application_variables.cpp

  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_cross_section.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/static_condensation_utility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/formfinding_io_utility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp

  ## RESPONSE UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/adjoint_structural_response_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/adjoint_local_stress_response_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/adjoint_nodal_displacement_response_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/adjoint_linear_strain_energy_response_function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/stress_response_definitions.cpp

  ## DIFFERENTIATION UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/response_utilities/element_finite_difference_utility.cpp

  ## PROCESSES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/postprocess_eigenvalues_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/total_structural_mass_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/prism_neighbours_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/shell_to_solid_shell_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/solid_shell_thickness_compute_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/spr_error_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/impose_rigid_movement_process.cpp

  ## RESPONSE FUNCTIONS
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_processes/replace_elements_and_conditions_for_adjoint_problem_process.cpp

  ##CONSTITUTIVE LAW
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/truss_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/truss_plasticity_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/beam_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_stress.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_strain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/elastic_isotropic_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/axisym_elastic_isotropic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_kirchhoff_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_stress_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_kirchhoff_plane_strain_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_orthotropic_2D_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_j2_plasticity_plane_strain_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_j2_plasticity_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_isotropic_damage_plane_strain_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_isotropic_damage_3D_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_isotropic_damage_factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_isotropic_plasticity_factory.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generic_small_strain_isotropic_damage.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generic_small_strain_isotropic_plasticity.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generic_small_strain_viscoplasticity_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/viscous_generalized_kelvin.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/viscous_generalized_maxwell.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generic_small_strain_d_plus_d_minus_damage.cpp

  ## ELEMENTS
  #Adding truss element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/truss_element_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/truss_element_linear_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cable_element_3D2N.cpp

  # Adding the beam element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cr_beam_element_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cr_beam_element_linear_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cr_beam_element_2D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cr_beam_element_linear_2D2N.cpp

  # Adding shells and membranes elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/isotropic_shell_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_shell_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_element_3D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thin_element_3D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thin_element_3D3N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_element_3D3N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/prestress_membrane_element.cpp

  # Adding elements for adjoint sensitivity analysis
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_elements/adjoint_finite_difference_base_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_elements/adjoint_finite_difference_shell_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_elements/adjoint_finite_difference_cr_beam_element_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_3D2N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_elements/adjoint_finite_difference_truss_element_linear_3D2N.cpp

  # Solid elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_solid_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_small_displacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_total_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_bbar.cpp

  # Adding the SPRISM element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid_shell_element_sprism_3D6N.cpp

  # Adding the nodal concentrated element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nodal_concentrated_element.cpp

  # Adding the spring + damper element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spring_damper_element_3D2N.cpp

  ## CONDITIONS
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/base_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_contact_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_point_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_load_condition_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_line_load_condition_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_condition_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_moment_condition_3d.cpp

  # Adding conditions for adjoint sensitivity analysis
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_conditions/adjoint_semi_analytic_point_load_condition.cpp
)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE
  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/structural_mechanics_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_response_functions_to_python.cpp
)

add_library(KratosStructuralMechanicsCore SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosStructuralMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStructuralMechanicsApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosStructuralMechanicsCore)
    cotire(KratosStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
