set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosContactMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication )

################### PYBIND11
include(pybind11Tools)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication)

if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle)
  include_directories( ${TRIANGLE_DIR} )
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication/external_libraries/tetgen)
  include_directories( ${TETGEN_DIR} )
  link_directories( ${TETGEN_DIR} )
endif(NOT DEFINED ${INCLUDE_TETGEN})

# Generate variables with the sources
set( KRATOS_CONTACT_MECHANICS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application_variables.cpp

	## deformable contact conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/contact_domain_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/contact_domain_LM_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/contact_domain_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/axisym_contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/deformable_contact/axisym_contact_domain_penalty_2D_condition.cpp

	## thermal contact conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_contact/thermal_contact_domain_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_contact/thermal_contact_domain_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_contact/axisym_thermal_contact_domain_penalty_2D_condition.cpp

	## rigid contact conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/point_rigid_contact_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/point_rigid_contact_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/axisym_point_rigid_contact_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/EP_point_rigid_contact_penalty_wP_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/EP_axisym_point_rigid_contact_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact/rigid_body_point_rigid_contact_condition.cpp

	## hydraulic contact conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/hydraulic_contact/hydraulic_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/hydraulic_contact/hydraulic_axisym_rigid_contact_penalty_2D_condition.cpp


	## beam contact conditions
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_conctact/beam_point_rigid_contact_condition.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_conctact/beam_point_rigid_contact_penalty_3D_condition.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_conctact/beam_point_rigid_contact_LM_3D_condition.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/beam_conctact/beam_point_pressure_condition.cpp

	## rigid body links
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_body_links/rigid_body_point_link_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_body_links/rigid_body_point_link_segregated_V_condition.cpp

	##elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rigid_body_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rigid_body_segregated_V_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/translatory_rigid_body_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/translatory_rigid_body_segregated_V_element.cpp

	##friction laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/coulomb_adhesion_friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/hardening_coulomb_friction_law.cpp


	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/contact_domain_utilities.cpp

	##meshers
	${CMAKE_CURRENT_SOURCE_DIR}/custom_meshers/contact_domain_3D_mesher.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_meshers/contact_domain_2D_mesher.cpp

)

## Generate variables with the sources
set( KRATOS_CONTACT_MECHANICS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_meshers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_friction_laws_to_python.cpp
)

add_library(KratosContactMechanicsCore SHARED ${KRATOS_CONTACT_MECHANICS_APPLICATION_CORE})
target_link_libraries(KratosContactMechanicsCore PUBLIC KratosCore KratosDelaunayMeshingCore)
set_target_properties(KratosContactMechanicsCore PROPERTIES COMPILE_DEFINITIONS "CONTACT_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosContactMechanicsApplication MODULE THIN_LTO ${KRATOS_CONTACT_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosContactMechanicsApplication PRIVATE KratosContactMechanicsCore)
set_target_properties(KratosContactMechanicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosContactMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosContactMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosContactMechanicsCore)
    cotire(KratosContactMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ContactMechanicsApplication.py KratosMultiphysics/ContactMechanicsApplication/__init__.py )

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
install(TARGETS KratosContactMechanicsCore DESTINATION libs )
install(TARGETS KratosContactMechanicsApplication DESTINATION libs )