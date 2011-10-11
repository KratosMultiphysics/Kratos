set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMetisApplication ****")
 
include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#find mkl solver package -- this is our own find package!!
find_package(PARMETIS REQUIRED)
include_directories(${PARMETIS_INCLUDE_DIR})
link_directories(${PARMETIS_LIBRARIES})

##add definitions needed
add_definitions(-DMKL_ILP64)

## generate variables with the sources
set( KRATOS_METIS_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/metis_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_metis_python_application.cpp 
)


# ###############################################################
add_library(KratosMetisApplication SHARED ${KRATOS_METIS_APPLICATION_SOURCES})
target_link_libraries(KratosMetisApplication KratosCore ${METIS_LIBRARY} )
set_target_properties(KratosMetisApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMetisApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosMetisApplication DESTINATION libs )

