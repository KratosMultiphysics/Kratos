set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMyLaplacianApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_MY_LAPLACIAN_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/my_laplacian_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/my_laplacian_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/my_laplacian_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/my_laplacian_element.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosMyLaplacianApplication SHARED ${KRATOS_MY_LAPLACIAN_APPLICATION_SOURCES})
target_link_libraries(KratosMyLaplacianApplication KratosCore )
set_target_properties(KratosMyLaplacianApplication PROPERTIES PREFIX "")
install(TARGETS KratosMyLaplacianApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMyLaplacianApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMyLaplacianApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MyLaplacianApplication.py" DESTINATION KratosMultiphysics )
