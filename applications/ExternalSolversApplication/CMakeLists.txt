set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExternalSolversApplication ****")
message("KratosExternalSolversApplication not yet ported to cmake")
# add_definitions( -DTETGEN143 )
# add_definitions( -DTETLIBRARY )
# 
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.1 )

# ##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.1)

## generate variables with the sources
set( KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/externalsolvers_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/externalsolvers_python_application.cpp 
)
# 
# ###############################################################
# ## define library Kratos which defines the basic python interface
add_library(KratosExternalSolversApplication SHARED ${KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES})
target_link_libraries(KratosExternalSolversApplication KratosCore triangle tetgen msuite )
set_target_properties(KratosExternalSolversApplication PROPERTIES PREFIX "")
install(TARGETS KratosExternalSolversApplication DESTINATION libs )

