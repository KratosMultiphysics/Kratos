set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSwimmingDEMApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)

include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Dependecies on other applications
include_directories( ${KRATOS_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/SwimmingDEMApplication )

## Generate variables with the sources
set( KRATOS_SWIMMING_DEM_APPLICATION_CORE_SOURCES
    ## Application class
    ${CMAKE_CURRENT_SOURCE_DIR}/swimming_DEM_application.cpp

    ## Strategies
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/hybrid_bashforth_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/terminal_velocity_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/symplectic_euler_old_velocity_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/strategies/adams_bashforth_strategy.cpp

    ## Constitutive Laws
    # Hydrodynamic Interaction Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hydrodynamic_interaction_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/power_law_hydrodynamic_interaction_law.cpp

    # Buoyancy laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/buoyancy_laws/buoyancy_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/buoyancy_laws/archimedes_buoyancy_law.cpp

    # Drag Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/stokes_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/beetstra_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/schiller_and_naumann_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/haider_and_levenspiel_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/chien_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/ganser_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/shah_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/newton_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/dallavalle_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/difelice_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/rong_drag_law.cpp

    # Virtual mass Force Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/virtual_mass_force_laws/virtual_mass_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/virtual_mass_force_laws/auton_hunt_prudhomme_virtual_mass_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/virtual_mass_force_laws/zuber_virtual_mass_force_law.cpp

    # Undisturbed Force Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/undisturbed_force_laws/undisturbed_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/undisturbed_force_laws/pressure_gradient_undisturbed_force_law.cpp

    # History Force Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/history_force_laws/history_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/history_force_laws/boussinesq_basset_history_force_law.cpp

    # Vorticity-Induced Lift Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/vorticity_induced_lift_laws/vorticity_induced_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/vorticity_induced_lift_laws/saffman_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/vorticity_induced_lift_laws/mei_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/vorticity_induced_lift_laws/el_samni_lift_law.cpp

    # Rotation-Induced Lift Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rotation_induced_lift_laws/rotation_induced_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rotation_induced_lift_laws/rubinow_and_keller_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rotation_induced_lift_laws/oesterle_dinh_lift_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rotation_induced_lift_laws/loth_rotation_induced_lift_law.cpp

    # Steady Viscous torque laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/steady_viscous_torque_laws/steady_viscous_torque_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/steady_viscous_torque_laws/rubinow_and_keller_torque_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/steady_viscous_torque_laws/loth_steady_viscous_torque_law.cpp

    #Process
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/bump_transient_porosity_solution_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/porosity_solution_and_sinusoidal_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/skrzypacz_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/skrzypacz_homogeneous_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/codina2001_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/sinusoidal_porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/sinusoidal_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/hyperbolic_tangential_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/hyperbolic_tangential_porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/plateau_bump_porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/plateau_bump_porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/plateau_bump_transient_porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/non_divergence_free_transient_porosity_solution_transient_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/porosity_solution_and_body_force_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/flow_past_cylinder_porosity_field_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/plateau_linear_porosity_field_process.cpp

    ## Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_dem_coupled.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_dem_coupled_weak.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_laplacian_simplex_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_mat_deriv_simplex_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_component_gradient_simplex_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_gradient_Pouliot_2012.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_gradient_Pouliot_2012_edge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_velocity_laplacian_component.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_velocity_laplacian.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/calculate_fluid_fraction_gradient_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_rigid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/swimming_particle.cpp

    ## Conditions
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_dem_coupled_wall_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/calculate_laplacian_simplex_condition.cpp

    ## Utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/binbased_DEM_fluid_coupled_mapping.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/basset_force_tools.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/velocity_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/cellular_flow_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/ethier_flow_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/product_of_sines_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/pouliot_flow_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/pouliot_flow_field_2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/constant_velocity_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/statistics/sampling_tool.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/sets/space_time_set.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/derivative_recovery.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/derivative_recovery_meshing_tools.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/field_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/inlets/bentonite_force_based_inlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/fluid_field_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fields/shear_flow_1D_with_exponential_viscosity_field.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/stationarity_check.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/error_norm_calculator_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/averaging_variables_utility.cpp
)

## Define python interface sources
set( KRATOS_SWIMMING_DEM_APPLICATION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/swimming_DEM_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosSwimmingDEMCore SHARED ${KRATOS_SWIMMING_DEM_APPLICATION_CORE_SOURCES})
target_link_libraries(KratosSwimmingDEMCore KratosDEMCore KratosFluidDynamicsCore KratosCore)
set_target_properties(KratosSwimmingDEMCore PROPERTIES COMPILE_DEFINITIONS "SWIMMING_DEM_APPLICATION=EXPORT,API")

## python module
pybind11_add_module(KratosSwimmingDEMApplication MODULE ${KRATOS_SWIMMING_DEM_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosSwimmingDEMApplication PUBLIC KratosSwimmingDEMCore)
set_target_properties(KratosSwimmingDEMApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosSwimmingDEMCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosSwimmingDEMApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/SwimmingDEMApplication.py KratosMultiphysics/SwimmingDEMApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/tests KratosMultiphysics/${CURRENT_DIR_NAME} )

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

## Install targets
install(TARGETS KratosSwimmingDEMCore DESTINATION libs )
install(TARGETS KratosSwimmingDEMApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosSwimmingDEMCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosSwimmingDEMApplication" PARENT_SCOPE)
