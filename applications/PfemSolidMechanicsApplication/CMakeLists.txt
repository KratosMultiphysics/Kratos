set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemSolidMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ContactMechanicsApplication )


if(NOT (${SOLID_MECHANICS_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Solid Mechanics Application must be compiled before this one. Set SOLID_MECHANICS_APPLICATION=ON in the configure file.")
endif(NOT (${SOLID_MECHANICS_APPLICATION} MATCHES "ON"))

if(NOT (${PFEM_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Pfem Base Application must be compiled before this one. Set PFEM_APPLICATION=ON in the configure file.")
endif(NOT (${PFEM_APPLICATION} MATCHES "ON"))

if(NOT (${CONTACT_MECHANICS_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Contact Mechanics Application must be compiled before this one. Set CONTACT_MECHANICS_APPLICATION=ON in the configure file.")
endif(NOT (${CONTACT_MECHANICS_APPLICATION} MATCHES "ON"))


if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
  include_directories( ${TRIANGLE_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../PfemApplication/external_modules/tetgen/)
  include_directories( ${TETGEN_INSTALL_DIR} )
  link_directories( ${TETGEN_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TETGEN})

## generate variables with the sources
set( KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp

	##elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_updated_lagrangian_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_updated_lagrangian_U_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_wP_Stab_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_wP_Stab_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_W_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_W_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_W_wP_DME_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_HO_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_FIC_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_DME_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_stab_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_U_W_wP_element.cpp

	## mixed elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_P_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_Pressure_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_P_wP_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_J_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_J_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_J_W_wP_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_J_W_wP_DME_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_Pressure_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian_U_P_wP_element.cpp


	##conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_rigid_contact_LM_3D_condition.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_point_pressure_condition.cpp


	##constitutive laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/non_linear_hencky_plastic_U_P_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_cam_clay_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_casm_axisym_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_casm_cem_axisym_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/borja_hencky_casm_cem_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.cpp


    ##hardening laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/casm_cem_hardening_law.cpp

	##yield criteria
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/casm_yield_criterion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/casm_cem_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/J2_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/tresca_yield_criterion.cpp

	##flow rules
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/borja_casm_explicit_plastic_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/borja_casm_cem_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.cpp


	##processes
	${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/set_mechanical_initial_state_process.cpp

	##utilities
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/contact_domain_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/water_pressure_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/water_pressure_utilities_Jacobian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/axisym_water_pressure_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/axisym_water_pressure_utilities_Jacobian.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/advanced_mesh_data_transfer_utilities.cpp


)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemSolidMechanicsApplication SHARED ${KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore)
target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore KratosSolidMechanicsApplication KratosContactMechanicsApplication)
set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES PREFIX "")
set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "PFEM_SOLID_MECHANICS_APPLICATION=EXPORT,API")

install(TARGETS KratosPfemSolidMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemSolidMechanicsApplication.py" DESTINATION KratosMultiphysics )
