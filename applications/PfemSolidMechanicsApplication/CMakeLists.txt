set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemSolidMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ContactMechanicsApplication )

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/SolidMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ContactMechanicsApplication)

if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
  include_directories( ${TRIANGLE_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication/external_modules/tetgen/)
  include_directories( ${TETGEN_INSTALL_DIR} )
  link_directories( ${TETGEN_INSTALL_DIR} )
endif(NOT DEFINED ${INCLUDE_TETGEN})

## Generate variables with the sources
set( KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application_variables.cpp

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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_J_W_wP_DME_element.cpp

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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/new_hencky_tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_tresca_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/new_hencky_tresca_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.cpp


    ##hardening laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.cpp

	##yield criteria
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/J2_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/new_tresca_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/tresca_yield_criterion.cpp

	##flow rules
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.cpp

	##processes
	${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/set_mechanical_initial_state_process.cpp

	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/water_pressure_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/water_pressure_utilities_Jacobian.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/axisym_water_pressure_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/axisym_water_pressure_utilities_Jacobian.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/advanced_mesh_data_transfer_utilities.cpp

)

## Generate variables with the sources
set( KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

add_library(KratosPfemSolidMechanicsCore SHARED ${KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_CORE})
target_link_libraries(KratosPfemSolidMechanicsCore PUBLIC KratosCore KratosSolidMechanicsCore KratosDelaunayMeshingCore KratosContactMechanicsCore)
#target_link_libraries(KratosPfemSolidMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosPfemSolidMechanicsCore PROPERTIES COMPILE_DEFINITIONS "PFEM_SOLID_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPfemSolidMechanicsApplication MODULE THIN_LTO ${KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPfemSolidMechanicsApplication PRIVATE KratosPfemSolidMechanicsCore)
set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosPfemSolidMechanicsCore)
    cotire(KratosPfemSolidMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PfemSolidMechanicsApplication.py KratosMultiphysics/PfemSolidMechanicsApplication/__init__.py )

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
install(TARGETS KratosPfemSolidMechanicsCore DESTINATION libs )
install(TARGETS KratosPfemSolidMechanicsApplication DESTINATION libs )