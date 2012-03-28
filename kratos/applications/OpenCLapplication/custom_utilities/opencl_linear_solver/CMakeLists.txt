set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring Test OpenCL solver ****")

#find mkl solver package -- this is our own find package!!
# find_package(OPENCL REQUIRED)
# include_directories(${OPENCL_INCLUDE_DIRS})
# link_directories(${OPENCL_LIBRARIES})

## generate variables with the sources
set( OPENCL_TEST_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/main.cpp 
)

###############################################################
## define library Kratos which defines the basic python interface
add_executable(test_opencl_solver ${OPENCL_TEST_SOURCES})
target_link_libraries(test_opencl_solver KratosCore ${OPENCL_LIBRARIES} rt )
install(TARGETS test_opencl_solver DESTINATION applications/OpenCLapplication/custom_utilities/opencl_linear_solver )