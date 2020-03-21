set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosParticleMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_PARTICLE_MECHANICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/particle_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/particle_mechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/particle_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_UP.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_quadrilateral.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_axisymmetry.cpp

	###constitutive laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_UP_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plane_strain_UP_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_UP_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_plane_strain_UP_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_UP_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_plane_strain_UP_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_strain_softening_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_strain_softening_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_borja_cam_clay_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_borja_cam_clay_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_borja_cam_clay_axisym_2D_law.cpp

	###yield criteria
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/yield_criteria/mc_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.cpp

	###yield criteria
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hardening_laws/exponential_strain_softening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hardening_laws/cam_clay_hardening_law.cpp

	###flow rules
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/particle_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/mc_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/mc_strain_softening_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.cpp

	###conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_base_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_axisym_point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_line_load_condition_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_axisym_line_load_condition_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/grid_based_conditions/mpm_grid_surface_load_condition_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_base_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_base_dirichlet_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_penalty_dirichlet_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_penalty_coupling_interface_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_base_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/particle_based_conditions/mpm_particle_point_load_condition.cpp

	###utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mpm_particle_generator_utility.cpp

)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE  KRATOS_PARTICLE_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
pybind11_add_module(KratosParticleMechanicsApplication MODULE THIN_LTO ${KRATOS_PARTICLE_MECHANICS_APPLICATION_SOURCES} ${KRATOS_PARTICLE_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosParticleMechanicsApplication PUBLIC KratosCore)
set_target_properties(KratosParticleMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "PARTICLE_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosParticleMechanicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosParticleMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosParticleMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosParticleMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ParticleMechanicsApplication.py KratosMultiphysics/ParticleMechanicsApplication/__init__.py )

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
install(TARGETS KratosParticleMechanicsApplication DESTINATION libs )