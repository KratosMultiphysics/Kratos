set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDEMApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set(KRATOS_DEM_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/DEM_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/strategies/explicit_solver_strategy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/strategies/explicit_solver_continuum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spheric_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nanoparticle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/analytic_spheric_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spheric_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ice_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_spheric_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/sintering_spheric_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/bonding_spheric_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/Particle_Contact_Element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cylinder_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cylinder_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cluster3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rigid_body_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ship_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/contact_info_spheric_particle.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/analytic_rigid_body_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/RigidFace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/analytic_RigidFace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/RigidEdge.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/SolidFace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dem_wall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mapping_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_discontinuum_constitutive_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_continuum_constitutive_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/singlespherecluster3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/dem_integration_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/forward_euler_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/symplectic_euler_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/taylor_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/velocity_verlet_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/runge_kutta_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/quaternion_integration_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/glued_to_wall_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_2D_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_torque_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_dev_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_2D_dev_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_fabric_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_soft_torque_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_soft_torque_with_noise_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_Rankine_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_Mohr_Coulomb_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_CamClay_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/dem_kdem_fissured_rock_cl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_ExponentialHC_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_HighStiffness_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_confined_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_confined_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_sintering_continuum_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Bentonite_Colloid_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_JKR_cohesive_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_DMT_cohesive_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Stress_dependent_cohesive_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/dem_d_linear_custom_constants_cl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Conical_damage_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/dem_kdem_2d_cl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/dem_kdem_fabric_2d_cl.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_viscous_Coulomb_Nestle_CL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/inlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/force_based_inlet.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/create_and_destroy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/properties_proxies.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/dem_fem_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/excavator_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/analytic_tools/analytic_model_part_filler.cpp
####    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/analytic_tools/analytic_particle_watcher.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/analytic_tools/particles_history_watcher.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_search/bins_dynamic_objects_periodic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/move_mesh_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/stationarity_checker.cpp
)

set(KRATOS_DEM_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/DEM_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/AuxiliaryUtilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/analytic_tools/analytic_face_watcher.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/analytic_tools/analytic_particle_watcher.cpp
)

add_library(KratosDEMCore SHARED ${KRATOS_DEM_APPLICATION_CORE})
target_link_libraries(KratosDEMCore PUBLIC KratosCore)
set_target_properties(KratosDEMCore PROPERTIES COMPILE_DEFINITIONS "DEM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDEMApplication MODULE THIN_LTO ${KRATOS_DEM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDEMApplication PUBLIC KratosDEMCore)
set_target_properties(KratosDEMApplication PROPERTIES PREFIX "")


if(${ACTIVATE_DEBUG_MACRO} MATCHES ON)        #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDEMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.clu" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosDEMApplication subdir inc_dirs = ${inc_dirs}")
# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosDEMCore)
  cotire(KratosDEMApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosDEMCore DESTINATION libs )
install(TARGETS KratosDEMApplication DESTINATION libs )

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DEMApplication.py" DESTINATION "KratosMultiphysics/DEMApplication" RENAME "__init__.py" )
