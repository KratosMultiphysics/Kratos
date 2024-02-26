set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPoromechanicsApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/LinearSolversApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/PoromechanicsApplication )

## Define application core sources (pure C++ parts)
set( KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES

    ## Application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application_variables.cpp

    # Utilities

    # Processes

    # Constitutive Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_flow_rules/local_damage_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_flow_rules/nonlocal_damage_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastoplastic_mohr_coulomb_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastoplastic_mod_mohr_coulomb_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/isotropic_damage_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/isotropic_damage_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastic_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/elastic_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/bilinear_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/bilinear_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/exponential_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/interface_element_laws/exponential_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/simo_ju_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/modified_mises_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/history_linear_elastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/history_linear_elastic_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/history_linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_hardening_laws/exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_hardening_laws/modified_exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_yield_criteria/simo_ju_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_yield_criteria/modified_mises_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_flow_rules/flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/custom_flow_rules/isotropic_damage_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/hyperelastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/linear_elastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/linear_elastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/hyperelastic_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/continuum_laws/linear_elastic_plastic_3D_law.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/U_Pl_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/U_Pl_small_strain_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/U_Pl_small_strain_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/U_Pl_small_strain_link_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/U_Pl_small_strain_FIC_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/one-phase_flow/small_strain_U_Pl_diff_order_element.cpp

    # Conditions
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_force_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_normal_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_liquid_discharge_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_FIC_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_face_load_interface_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/U_Pl_normal_liquid_flux_interface_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/general_U_Pl_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/line_load_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/line_normal_load_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/line_normal_liquid_flux_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/surface_load_3D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/surface_normal_load_3D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/one-phase_flow/surface_normal_liquid_flux_3D_diff_order_condition.cpp
)

## Define python interface sources
set( KRATOS_POROMECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/poromechanics_python_application.cpp
)

## PoromechanicsApplication core library (C++ parts)
add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
target_link_libraries( KratosPoromechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## PoromechanicsApplication core library (C++ parts)
#add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
#target_link_libraries( KratosPoromechanicsCore PRIVATE KratosCore KratosStructuralMechanicsCore)
set_target_properties( KratosPoromechanicsCore PROPERTIES COMPILE_DEFINITIONS "POROMECHANICS_APPLICATION=EXPORT,API")

## PoromechanicsApplication python module
pybind11_add_module(KratosPoromechanicsApplication MODULE THIN_LTO ${KRATOS_POROMECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosPoromechanicsApplication PRIVATE KratosPoromechanicsCore)
set_target_properties(KratosPoromechanicsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosPoromechanicsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosPoromechanicsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PoromechanicsApplication.py KratosMultiphysics/PoromechanicsApplication/__init__.py )

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

# Install targets
install(TARGETS KratosPoromechanicsCore DESTINATION libs )
install(TARGETS KratosPoromechanicsApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosPoromechanicsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosPoromechanicsApplication" PARENT_SCOPE)