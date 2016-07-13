set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosContactMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemBaseApplication )

# generate variables with the sources
set( KRATOS_CONTACT_MECHANICS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/contact_mechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

	##conditions
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_condition.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rigid_contact_condition.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosContactMechanicsApplication SHARED ${KRATOS_CONTACT_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosContactMechanicsApplication KratosCore KratosPfemBaseApplication)
set_target_properties(KratosContactMechanicsApplication PROPERTIES PREFIX "")
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
