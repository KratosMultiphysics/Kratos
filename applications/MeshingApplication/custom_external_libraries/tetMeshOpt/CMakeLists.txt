set(CMAKE_INCLUDE_CURRENT_DIR ON)

##removing optimizations

## generate variables with the sources
set( MESH_OPT      	
    ${CMAKE_CURRENT_SOURCE_DIR}/Math3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_delphiClasses.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_elementCluster.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_MeshLoaders.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ParallelFunctions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ProcessTime.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_qualityMetrics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ShowMetrics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_TetraFunctions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_tools.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/u_Types.cpp	
     )

###############################################################
##removing optimizations

add_library(tetMeshOpt STATIC ${MESH_OPT}  )
target_link_libraries(tetMeshOpt  )
install(TARGETS tetMeshOpt DESTINATION libs )

###############################################################
