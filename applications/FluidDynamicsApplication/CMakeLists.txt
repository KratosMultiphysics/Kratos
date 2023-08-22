set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
set( KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES

  # application class and variables
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application_variables.cpp

  # utilities (compiled first because they are used within the element cpps)
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/estimate_dt_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_auxiliary_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_calculation_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_characteristic_numbers_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_element_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_element_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/drag_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/statistics_record.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/vorticity_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/acceleration_limitation_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_adjoint_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_test_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_adjoint_test_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_adjoint_slip_utilities.cpp

  # elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_adjoint_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qs_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qs_vms_dem_coupled.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/alternative_qs_vms_dem_coupled.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/d_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/d_vms_dem_coupled.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/alternative_d_vms_dem_coupled.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/symbolic_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_fluid_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_fluid_element_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dynamic_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spalart_allmaras.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stationary_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stokes_3D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stokes_3D_twofluid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/navier_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_ausas_navier_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes_explicit_2D3N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes_explicit_2D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes_explicit_3D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_fluid_navier_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_fluid_navier_stokes_alpha_method.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms_adjoint_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/distance_smoothing_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/data_containers/qs_vms/qs_vms_derivative_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/data_containers/qs_vms/qs_vms_residual_derivatives.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/weakly_compressible_navier_stokes.cpp

  # conditions
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_werner_wengle_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_generalized_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_periodic_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/navier_stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/two_fluid_navier_stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/embedded_ausas_navier_stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/adjoint_monolithic_wall_condition.cpp

  # fluid constitutive laws
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/fluid_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bingham_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/euler_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/euler_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_temperature_dependent_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_temperature_dependent_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_two_fluid_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_two_fluid_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/herschel_bulkley_3d_law.cpp

  # processes
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/apply_compressible_navier_stokes_boundary_conditions_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/boussinesq_force_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/calulate_levelset_consistent_nodal_gradient_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_pressure_coefficient_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/distance_modification_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/distance_smoothing_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/embedded_skin_visualization_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mass_conservation_check_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/two_fluids_inlet_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/shock_capturing_entropy_viscosity_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/shock_capturing_physics_based_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_y_plus_process.cpp
)

## define python interface sources
set( KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_response_functions_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_dynamics_python_application.cpp
)

## define fluid dynamics test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FluidDynamicsApplication core library (C++ parts)
add_library( KratosFluidDynamicsCore SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES} )
target_link_libraries( KratosFluidDynamicsCore PUBLIC KratosCore)
set_target_properties( KratosFluidDynamicsCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_APPLICATION=EXPORT,API")

## FluidDynamicsApplication python module
pybind11_add_module( KratosFluidDynamicsApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsApplication PRIVATE KratosFluidDynamicsCore)
set_target_properties( KratosFluidDynamicsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidDynamicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFluidDynamicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsApplication.py KratosMultiphysics/FluidDynamicsApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/FluidDynamicsApplication/python_registry_lists.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
      )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosFluidDynamicsCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )

## If trilinos is avaliable, switch on the Fluid trilinos extension
if((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))
  message(STATUS "Enabling trilinos extension for KratosFluidDynamicsApplication")
  add_subdirectory(trilinos_extension)
endif((${USE_MPI} MATCHES ON) AND (${TRILINOS_FOUND}))

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidDynamicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidDynamicsApplication" PARENT_SCOPE)
