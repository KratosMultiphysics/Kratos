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
install(TARGETS KratosALEApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")