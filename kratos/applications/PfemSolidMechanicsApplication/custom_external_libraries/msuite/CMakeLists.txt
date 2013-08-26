set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_definitions( -D_NO_DELETE_NODES )

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src/util )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src/edt )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/src/malla )

## generate variables with the sources
set( MSUITE_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/case.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/cone.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/filename.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/octree.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/pline.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/punto.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/util/tiempo.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/alpha.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/edt.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/esfera.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/hydir.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/s_squeeze.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/slivers.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/vor_out.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/edt/voronoi.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/corte.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/delaunay.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/elemento.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/fforma.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/grabamalla.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/leemalla.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/malla.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/meshelm.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/no_remalla.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/orienta.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/refina.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/smooth.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/surface.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/src/malla/vecinos_y_frontera.cpp 
    )

###############################################################
add_library(msuite STATIC ${MSUITE_SOURCES}  )
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(msuite m)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# install(TARGETS msuite DESTINATION libs )
