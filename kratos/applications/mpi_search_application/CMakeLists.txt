set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMPISerchApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_MPI_SEARCH_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/mpi_search_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_mpi_search_python_application.cpp 
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMPISearchApplication SHARED ${KRATOS_MPI_SEARCH_APPLICATION_SOURCES})
target_link_libraries(KratosMPISearchApplication KratosCore )
set_target_properties(KratosMPISearchApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMPISearchApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install(TARGETS KratosMPISearchApplication DESTINATION libs )


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MPISearchApplication.py" DESTINATION KratosMultiphysics )
