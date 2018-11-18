set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMappingApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources

set( KRATOS_MAPPING_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/mapping_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mapping_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mappers_to_python.cpp
    # Mappers
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_mappers/nearest_neighbor_mapper.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_mappers/nearest_element_mapper.cpp
    # # Searching
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_searching/interface_communicator.cpp
    # # Utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/matrix_based_mapping_operation_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_preprocessor.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_factory.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_flags.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapper_utilities.cpp
)

set(KRATOS_MAPPING_APPLICATION_SOURCES_MPI "")
set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS "")
## If we have mpi support, use it
if(${MPI_NEEDED} MATCHES ON )
    include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application )
    find_package(TRILINOS REQUIRED)
    include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR}) # Included as "SYSTEM" to suppress warnings from it

    set( KRATOS_MAPPING_APPLICATION_SOURCES_MPI ${KRATOS_MAPPING_APPLICATION_SOURCES_MPI}
        # Searching
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_searching/interface_communicator_mpi.cpp
        # Utilities
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/matrix_based_mapping_operation_utility_mpi.cpp
    )

    set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES})
    set(KRATOS_MAPPING_APPLICATION_EXTRA_LIBS ${KRATOS_MAPPING_APPLICATION_EXTRA_LIBS} ${TRILINOS_LIBRARIES})
endif(${MPI_NEEDED} MATCHES ON )

## define test sources (cpp tests)
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_MAPPING_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMappingApplication MODULE THIN_LTO ${KRATOS_MAPPING_APPLICATION_SOURCES} ${KRATOS_MAPPING_APPLICATION_SOURCES_MPI} ${KRATOS_MAPPING_TESTING_SOURCES})
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
