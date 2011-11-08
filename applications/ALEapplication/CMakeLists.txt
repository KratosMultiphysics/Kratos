set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosALEApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_ALE_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/laplacian_meshmoving_element_2d.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/laplacian_meshmoving_element_3d.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/ale_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/ale_python_application.cpp 
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosALEApplication SHARED ${KRATOS_ALE_APPLICATION_SOURCES})
target_link_libraries(KratosALEApplication KratosCore )
set_target_properties(KratosALEApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosALEApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosALEApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")
