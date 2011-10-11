set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMKLSolversApplication ****")
 
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/external_includes )

#find mkl solver package -- this is our own find package!!
find_package(MKLSOLVER REQUIRED)
include_directories(${MKLSOLVER_INCLUDE_DIR})
link_directories(${MKLSOLVER_LIBRARIES})
message("mkl include dir ${MKLSOLVER_INCLUDE_DIR}")
message("mkl libs ${MKLSOLVER_LIBRARIES}")

##add definitions needed
add_definitions(-DMKL_ILP64)
		  
## generate variables with the sources
set( KRATOS_MKL_SOLVERS_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/mkl_solvers_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mkl_solvers_python_application.cpp 
)


# ###############################################################
add_library(KratosMKLSolversApplication SHARED ${KRATOS_MKL_SOLVERS_APPLICATION_SOURCES})
target_link_libraries(KratosMKLSolversApplication KratosCore ${MKLSOLVER_LIBRARIES} )
set_target_properties(KratosMKLSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMKLSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosMKLSolversApplication DESTINATION libs )

