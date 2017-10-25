set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidDynamicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FLUID_DYNAMICS_APPLICATION_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_dynamics_application_variables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_dynamics_python_application.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dynamic_vms.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spalart_allmaras.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stationary_stokes.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stokes_3D.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/stokes_3D_twofluid.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/navier_stokes.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_werner_wengle_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_generalized_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition_discontinuous.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/stokes_wall_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fs_periodic_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/navier_stokes_wall_condition.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bingham_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_2d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/newtonian_3d_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/herschel_bulkey_3d_law.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/boussinesq_force_process.cpp

)

## fluid dynamics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FLUID_DYNAMICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFluidDynamicsApplication SHARED ${KRATOS_FLUID_DYNAMICS_APPLICATION_SOURCES} ${KRATOS_FLUID_DYNAMICS_TESTING_SOURCES})
target_link_libraries(KratosFluidDynamicsApplication KratosCore )
set_target_properties(KratosFluidDynamicsApplication PROPERTIES COMPILE_DEFINITIONS "FLUID_DYNAMICS_APPLICATION=EXPORT,API")
set_target_properties(KratosFluidDynamicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosFluidDynamicsApplication DESTINATION libs )

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FluidDynamicsApplication.py" DESTINATION KratosMultiphysics )
