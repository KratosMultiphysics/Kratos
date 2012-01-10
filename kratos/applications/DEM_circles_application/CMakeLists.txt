set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDEMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_DEM_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/DEM_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/2d_circular_particle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/spheric_particle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/DEM_python_application.cpp 
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosDEMApplication SHARED ${KRATOS_DEM_APPLICATION_SOURCES})
target_link_libraries(KratosDEMApplication KratosCore )
set_target_properties(KratosDEMApplication PROPERTIES PREFIX "")
install(TARGETS KratosDEMApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosDEMApplication subdir inc_dirs = ${inc_dirs}")