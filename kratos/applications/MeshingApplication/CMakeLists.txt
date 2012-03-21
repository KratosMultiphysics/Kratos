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
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)

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
target_link_libraries(KratosMeshingApplication KratosCore triangle tetgen msuite tetMeshOpt  )
set_target_properties(KratosMeshingApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMeshingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosMeshingApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MeshingApplication.py" DESTINATION KratosMultiphysics )
