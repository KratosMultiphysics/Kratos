set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRadiationApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_RADIATION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//radiation_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/radiation_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/radiation_3d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/radiation_2d.cpp


	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rad_face3D.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/rad_face2D.cpp




	

	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/your_new_element_2d.h
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosRadiationApplication MODULE THIN_LTO ${KRATOS_RADIATION_APPLICATION_SOURCES})
#add_library(KratosRadiationApplication SHARED ${KRATOS_RADIATION_APPLICATION_SOURCES})
target_link_libraries(KratosRadiationApplication KratosCore )
set_target_properties(KratosRadiationApplication PROPERTIES PREFIX "")
install(TARGETS KratosRadiationApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosRadiationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/RadiationApplication.py" DESTINATION KratosMultiphysics )
