set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPFEMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_PFEM_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/PFEM_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_pfem_python_application.cpp 
  	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_meshers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp      

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/free_surface_cond2d.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPFEMApplication SHARED ${KRATOS_PFEM_APPLICATION_SOURCES})
target_link_libraries(KratosPFEMApplication KratosCore )
set_target_properties(KratosPFEMApplication PROPERTIES PREFIX "")
install(TARGETS KratosPFEMApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")