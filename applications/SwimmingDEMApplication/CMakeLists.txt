set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSwimmingDEMApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)

include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## Dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SwimmingDEMApplication )

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
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/ganser_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/shah_drag_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/drag_laws/newton_drag_law.cpp

    # Inviscid Force Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/inviscid_force_laws/inviscid_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/inviscid_force_laws/auton_hunt_prudhomme_inviscid_force_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/inviscid_force_laws/zuber_inviscid_force_law.cpp

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

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosSwimmingDEMCore)
  cotire(KratosSwimmingDEMApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/SwimmingDEMApplication.py KratosMultiphysics/SwimmingDEMApplication/__init__.py )

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

## Install targets
install(TARGETS KratosSwimmingDEMCore DESTINATION libs )
install(TARGETS KratosSwimmingDEMApplication DESTINATION libs )