set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshingApplication ****")

add_definitions( -DTETGEN143 )
add_definitions( -DTETLIBRARY )

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )

##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)

## generate variables with the sources
set( KRATOS_MESHING_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/meshing_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_meshing_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_meshers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 	

	${CMAKE_CURRENT_SOURCE_DIR}/external_includes/mesh_suites.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMeshingApplication SHARED ${KRATOS_MESHING_APPLICATION_SOURCES})
target_link_libraries(KratosMeshingApplication KratosCore triangle tetgen msuite  )
set_target_properties(KratosMeshingApplication PROPERTIES PREFIX "")
install(TARGETS KratosMeshingApplication DESTINATION libs )

