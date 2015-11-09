set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosEmpireApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other application -- TRILINOS specific
include_directories( ${CMAKE_SOURCE_DIR}/applications/structural_application )

## generate variables with the sources
set( KRATOS_EMPIRE_APPLICATION_SOURCES     
        ${CMAKE_CURRENT_SOURCE_DIR}/empire_application.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/empire_python_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_SOURCE_DIR}/applications/structural_application/custom_elements/membrane_element.cpp

	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosEmpireApplication SHARED ${KRATOS_EMPIRE_APPLICATION_SOURCES})
target_link_libraries(KratosEmpireApplication KratosCore )
set_target_properties(KratosEmpireApplication PROPERTIES PREFIX "")
install(TARGETS KratosEmpireApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosEmpireApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosEmpireApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/EmpireApplication.py" DESTINATION KratosMultiphysics )
