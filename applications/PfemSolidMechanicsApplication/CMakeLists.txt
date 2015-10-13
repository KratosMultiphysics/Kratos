set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemSolidMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication/custom_external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication/custom_external_libraries/tetgen/tetgen1.4.3 )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)


## generate variables with the sources
set( KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/contact_domain_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/modeler_utilities.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/advanced_mesh_data_transfer_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_data_transfer_utilities.cpp
	
        
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/J2_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/tresca_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/mohr_coulomb_yield_criterion.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/matsuoka_nakai_yield_criterion.cpp


	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/cam_clay_explicit_plastic_flow_rule.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/linear_cam_clay_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/mohr_coulomb_explicit_plastic_flow_rule.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_associative_plastic_flow_rule.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/matsuoka_nakai_flow_rule.cpp
        

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_plane_strain_2D_law.cpp


	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_cam_clay_plane_strain_2D_law.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_cam_clay_axisym_2D_law.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_hencky_cam_clay_plane_strain_2D_law.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_hencky_cam_clay_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mohr_coulomb_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_mohr_coulomb_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_updated_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_updated_lagrangian_U_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_Stab_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_FIC_element.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_weight_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_wP_Stab_element.cpp
	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/custom_friction_laws/friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/custom_friction_laws/hardening_coulomb_friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/custom_friction_laws/coulomb_adhesion_friction_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/composite_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_contact_domain_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_point_rigid_contact_penalty_water_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_body_point_rigid_contact_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_LM_3D_condition.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_pressure_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/mesh_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/triangular_mesh_2D_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/contact_domain_2D_modeler.cpp


	${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/set_mechanical_initial_state_process.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemSolidMechanicsApplication SHARED ${KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore)
target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore KratosSolidMechanicsApplication KratosMeshingApplication triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore KratosSolidMechanicsApplication triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore triangle tetgen msuite tetMeshOpt)
set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosPfemSolidMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemSolidMechanicsApplication.py" DESTINATION KratosMultiphysics )
