set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFSIApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FSI_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/fsi_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fsi_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mappers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/AdvancedNMPointsMapper.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFSIApplication SHARED ${KRATOS_FSI_APPLICATION_SOURCES})
target_link_libraries(KratosFSIApplication KratosCore )
set_target_properties(KratosFSIApplication PROPERTIES PREFIX "")
install(TARGETS KratosFSIApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")