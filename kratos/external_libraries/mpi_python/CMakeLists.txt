set(CMAKE_INCLUDE_CURRENT_DIR ON)

## generate variables with the sources
set( MPI_PYTHON_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/python_interface.cpp
     )

add_definitions( -w )

###############################################################
add_library(mpipython SHARED ${MPI_PYTHON_SOURCES})

target_link_libraries(mpipython ${Boost_PYTHON_LIBRARIES} ${PYTHON_LIBRARIES} ${MPI_LIBRARIES})

set_target_properties(mpipython PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(mpipython PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# Install library and Python file
install(TARGETS mpipython DESTINATION libs )
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mpi.py" DESTINATION KratosMultiphysics )
