set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosOpenCLApplication ****")

#find mkl solver package -- this is our own find package!!
find_package(OPENCL REQUIRED)
include_directories(${OPENCL_INCLUDE_DIRS})
link_directories(${OPENCL_LIBRARIES})

add_definitions( -DVIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK )

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )

# message("OPENCL_INCLUDE_DIRS = ${OPENCL_INCLUDE_DIRS}")
# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")

## generate variables with the sources
set( KRATOS_OPENCL_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/OpenCLapplication.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/opencl_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/pugixml/src/pugixml.cpp 	
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosOpenCLApplication SHARED ${KRATOS_OPENCL_APPLICATION_SOURCES})
target_link_libraries(KratosOpenCLApplication KratosCore ${OPENCL_LIBRARIES} rt )
set_target_properties(KratosOpenCLApplication PROPERTIES PREFIX "")
install(TARGETS KratosOpenCLApplication DESTINATION libs )

