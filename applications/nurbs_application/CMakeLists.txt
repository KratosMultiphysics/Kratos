set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNurbsTestcaseApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
# include_directories( ${CMAKE_SOURCE_DIR}/kratos /home/matthias/Downloads/opennurbs5.0/ )

## generate variables with the sources
set( KRATOS_NURBSTESTCASE_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//nurbs_testcase_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/nurbs_testcase_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nurbs_poisson_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/NurbsModeler.cpp			             
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosNurbsTestcaseApplication SHARED ${KRATOS_NURBSTESTCASE_APPLICATION_SOURCES})
# set(OPENNURBS_LIB "/home/matthias/Downloads/opennurbs5.0/libopenNURBS.a")
target_link_libraries(KratosNurbsTestcaseApplication KratosCore ) #${OPENNURBS_LIB})
set_target_properties(KratosNurbsTestcaseApplication PROPERTIES PREFIX "")
install(TARGETS KratosNurbsTestcaseApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosNurbsTestcaseApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosNurbsTestcaseApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/NurbsTestcaseApplication.py" DESTINATION KratosMultiphysics )
