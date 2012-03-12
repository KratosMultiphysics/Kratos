set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_STRUCTURAL_APPLICATION_SOURCES  
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rigid_body_3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_anisotropic.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/beam_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_anisotropic_linear.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_isotropic.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/membrane_element.cpp 	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/linear_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/Hypoelastic_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/mixed_lagrangian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/kinematic_linear.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_2phase.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_1phase_small_strain.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_2phase_small_strain.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_3phase.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unsaturated_soils_element_3phase_small_strain.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/crisfield_truss_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ebst.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ebst_vel.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/pfem_contact_element3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/pfem_contact_element3D_vel.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/linear_elastic_truss.cpp
	
	${CMAKE_CURRENT_SOURCE_DIR}/structural_application.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_conditions_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/structural_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_parallel_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_fluency_criteria_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/hooks_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/external_isotropic_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_damage_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_damage_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_planestress_wrinkling.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_3d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/hyperelastic_3d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/hyperelastic_2d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/viscoelastic_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/viscofibers_2d.cpp # new!!
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/viscofibers_hypermatrix_2d.cpp # new!!
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/von_mises_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_2d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plane_strain.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/hypoelastic_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/fluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/external_isotropic_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/hooks_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/drucker_prager.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/casm.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/orthotropic_3d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plasticity_2d.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plane_stress_J2.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/mohr_coulomb_plane_strain.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/brittle_material_2d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/isotropic_rankine_damage_2d.cpp

	#rankine_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/isotropic_rankine_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/morh_coulomb_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/standard_morh_coulomb_yield_function.cpp
        #tresca_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/von_mises_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/modified_morh_coulomb_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/drucker_prager_yield_function.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/soft_hard_behavior/exponencial_softening.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/soft_hard_behavior/linear_softening.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/soft_hard_behavior/cohesion_softening.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/soft_hard_behavior/friction_softening.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/soft_hard_behavior/dilatancy_softening.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/fluency_criteria/energy_yield_function.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/master_contact_face_2d.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/master_contact_point_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/slave_contact_point_2d.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_segment_contact_link.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_point_contact_link.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointforce2D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointmoment3D.cpp	
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_vel_3D.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointforce3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face2D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face3D.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_pressure3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/faceforce3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/node_tying_lagrange.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/node_tying_lagrange_z.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/master_contact_face_3D_newmark.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/slave_contact_face_3D_newmark.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/master_contact_face_3D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/slave_contact_face_3D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_link_3D.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_link_3D_kinematic_linear.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_link_3D_newmark.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/slave_contact_point_3d.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_link_3D_explicit.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosStructuralApplication SHARED ${KRATOS_STRUCTURAL_APPLICATION_SOURCES})
target_link_libraries(KratosStructuralApplication KratosCore )
set_target_properties(KratosStructuralApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStructuralApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosStructuralApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralApplication.py" DESTINATION KratosMultiphysics )
