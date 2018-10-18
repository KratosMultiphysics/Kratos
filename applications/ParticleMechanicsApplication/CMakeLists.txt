set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosParticleMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )

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

	###constitutive laws
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_3d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_plane_strain_2d_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_UP_3d_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_plastic_plane_strain_UP_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_plane_strain_2D_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_plane_strain_UP_2D_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_3D_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_UP_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_strain_softening_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mc_strain_softening_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_borja_cam_clay_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_borja_cam_clay_plane_strain_2D_law.cpp

	###yield criteria
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/yield_criteria/mc_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/yield_criteria/modified_cam_clay_yield_criterion.cpp

	###yield criteria
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hardening_laws/exponential_strain_softening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hardening_laws/cam_clay_hardening_law.cpp

        ###flow rules
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/MPM_flow_rule.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/mc_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/mc_strain_softening_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/flow_rules/borja_cam_clay_plastic_flow_rule.cpp

	###conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mpm_base_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mpm_point_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mpm_line_load_condition_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mpm_surface_load_condition_3d.cpp
)

###############################################################
pybind11_add_module(KratosParticleMechanicsApplication MODULE THIN_LTO ${KRATOS_PARTICLE_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosParticleMechanicsApplication PUBLIC KratosCore KratosSolidMechanicsCore)
set_target_properties(KratosParticleMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "PARTICLE_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosParticleMechanicsApplication PROPERTIES PREFIX "")


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosParticleMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(USE_COTIRE MATCHES ON)
    cotire(KratosParticleMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosParticleMechanicsApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ParticleMechanicsApplication.py" DESTINATION KratosMultiphysics )
