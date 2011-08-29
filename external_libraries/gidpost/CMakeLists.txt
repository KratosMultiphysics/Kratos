set(CMAKE_INCLUDE_CURRENT_DIR ON)

## generate variables with the sources
set( GIDPOST_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/source/gidpost.cpp
 	${CMAKE_CURRENT_SOURCE_DIR}/source/gidpostInt.cpp
 	${CMAKE_CURRENT_SOURCE_DIR}/source/zlibint1.c
 	${CMAKE_CURRENT_SOURCE_DIR}/source/zlibint2.c
     )

add_definitions( -DUSE_CONST )

###############################################################
add_library(gidpost STATIC ${GIDPOST_SOURCES})

# install(TARGETS gidpost DESTINATION libs )