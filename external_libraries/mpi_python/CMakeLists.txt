set(CMAKE_INCLUDE_CURRENT_DIR ON)

include(pybind11Tools)

## generate variables with the sources
set( MPI_PYTHON_SOURCES
 	${CMAKE_CURRENT_SOURCE_DIR}/python_interface.cpp
     )

add_definitions( -w )

###############################################################
pybind11_add_module(mpipython MODULE THIN_LTO ${MPI_PYTHON_SOURCES})

target_link_libraries(mpipython PRIVATE ${MPI_LIBRARIES})

set_target_properties(mpipython PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(mpipython PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(mpipython PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(USE_COTIRE MATCHES ON)
    cotire(mpipython)
endif(USE_COTIRE MATCHES ON)
# Install library and Python file
install(TARGETS mpipython DESTINATION libs )
# install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mpi.py" DESTINATION KratosMultiphysics )
