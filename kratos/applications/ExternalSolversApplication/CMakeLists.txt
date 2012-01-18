set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExternalSolversApplication ****")
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
target_link_libraries(KratosExternalSolversApplication KratosCore superlu )
set_target_properties(KratosExternalSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

install(TARGETS KratosExternalSolversApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ExternalSolversApplication.py" DESTINATION KratosMultiphysics )
