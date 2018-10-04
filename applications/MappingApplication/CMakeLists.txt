set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMappingApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
# include_directories( ${CMAKE_SOURCE_DIR}/applications/mpi_search_application ) # TODO can be removed for now


# generate variables with the sources
set( KRATOS_MAPPING_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/mapping_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mapping_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mappers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_flags.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_factory.cpp
)

set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS "")
## If we have mpi support, use it
if(${MPI_NEEDED} MATCHES ON )
	set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES})
endif(${MPI_NEEDED} MATCHES ON )

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMappingApplication MODULE THIN_LTO ${KRATOS_MAPPING_APPLICATION_SOURCES})
target_link_libraries(KratosMappingApplication PRIVATE KratosCore ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosMappingApplication PROPERTIES PREFIX "")
install(TARGETS KratosMappingApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# enabling Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosMappingApplication)
endif(USE_COTIRE MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MappingApplication.py" DESTINATION KratosMultiphysics )
