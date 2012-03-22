set(CMAKE_INCLUDE_CURRENT_DIR ON)

##removing optimizations
SET_SOURCE_FILES_PROPERTIES(${UNOPTIMIZED_PREDICATES} PROPERTIES COMPILE_FLAGS -O0)
SET_SOURCE_FILES_PROPERTIES(${TETGEN_SOURCES} PROPERTIES COMPILE_FLAGS -O2)

## generate variables with the sources
set( MESH_OPT     
 	${CMAKE_CURRENT_SOURCE_DIR}/apiPrimitives.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Math3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_delphiClasses.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_elementCluster.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_MeshLoaders.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ParallelFunctions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ProcessTime.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_qualityMetrics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_ShowMetrics.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_TetraFunctions.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/u_Tools.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/u_Types.cpp	
     )

###############################################################
##removing optimizations
SET_SOURCE_FILES_PROPERTIES(${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS -O0)
add_definitions( -DTETOPTLIBRARY )

add_definitions( -w )

add_library(tetMeshOpt STATIC ${MESH_OPT}  )
target_link_libraries(tetMeshOpt  )
install(TARGETS tetMeshOpt DESTINATION libs )

###############################################################
