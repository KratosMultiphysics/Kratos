set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
set( KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES

  # application class and variables
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application_variables.cpp

  # utilities (compiled first because they are used within the element cpps)
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_element_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/fluid_element_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/drag_utilities.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/element_size_calculator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/vorticity_utilities.cpp

  # elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fluid_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qs_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/d_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/symbolic_navier_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_fluid_element.cpp
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
  #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes2D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_navier_stokes3D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/two_fluid_navier_stokes.cpp

  # conditions
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_werner_wengle_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_generalized_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_periodic_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/navier_stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/embedded_ausas_navier_stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/behr_wall_condition.cpp

  # fluid constitutive laws
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/fluid_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bingham_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/euler_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/euler_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_two_fluid_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_two_fluid_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/herschel_bulkley_3d_law.cpp

  # processes
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/boussinesq_force_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/distance_modification_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/embedded_skin_visualization_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/mass_conservation_check_process.cpp
)

## define python interface sources
set( KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_dynamics_python_application.cpp
)

## define fluid dynamics test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

#get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## FluidDynamicsApplication core library (C++ parts)
add_library( KratosFluidDynamicsCore SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_CORE_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES} )
target_link_libraries( KratosFluidDynamicsCore PUBLIC KratosCore)
set_target_properties( KratosFluidDynamicsCore PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_APPLICATION=EXPORT,API")

## FluidDynamicsApplication python module
pybind11_add_module( KratosFluidDynamicsApplication MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsApplication PRIVATE KratosFluidDynamicsCore)
set_target_properties( KratosFluidDynamicsApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosFluidDynamicsCore)
    cotire(KratosFluidDynamicsApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosFluidDynamicsCore DESTINATION libs )
install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsApplication.py" DESTINATION KratosMultiphysics )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

