set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosContactMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemApplication )


if(NOT (${PFEM_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The Pfem Base Application must be compiled before this one. Set PFEM_APPLICATION=ON in the configure file.")
endif(NOT (${PFEM_APPLICATION} MATCHES "ON"))

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


# generate variables with the sources
set( KRATOS_CONTACT_MECHANICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_friction_laws_to_python.cpp

	##conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_LM_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_contact_domain_LM_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_contact_domain_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_penalty_3D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_rigid_contact_penalty_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_point_rigid_contact_penalty_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_body_point_rigid_contact_condition.cpp

	##elements
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rigid_body_element.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/translatory_rigid_body_element.cpp

	##friction laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/hardening_coulomb_friction_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_friction/coulomb_adhesion_friction_law.cpp

	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/contact_domain_utilities.cpp

	##modelers
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/contact_domain_3D_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/contact_domain_2D_modeler.cpp

)

# define library Kratos which defines the basic python interface
add_library(KratosContactMechanicsApplication SHARED ${KRATOS_CONTACT_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosContactMechanicsApplication KratosCore KratosPfemApplication)
set_target_properties(KratosContactMechanicsApplication PROPERTIES PREFIX "")
set_target_properties(KratosContactMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "CONTACT_MECHANICS_APPLICATION=EXPORT,API")
install(TARGETS KratosContactMechanicsApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosContactMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosContactMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ContactMechanicsApplication.py" DESTINATION KratosMultiphysics )
