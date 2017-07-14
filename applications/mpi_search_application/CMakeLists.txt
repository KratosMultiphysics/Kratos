set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMPISerchApplication ****")

if(NOT (${DEM_APPLICATION} STREQUAL "ON"))
    message(SEND_ERROR "The DEM Application must be compiled before this one. Set DEM_APPLICATION=ON in the configure file.")
endif(NOT (${DEM_APPLICATION} STREQUAL "ON"))

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEM_application )
include_directories( ${CMAKE_SOURCE_DIR}/applications/mpi_search_application )

## generate variables with the sources
set( KRATOS_MPI_SEARCH_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/mpi_search_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_mpi_search_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

## Add mpi tests
if(${KRATOS_BUILD_TESTING} MATCHES ON)
	file(GLOB_RECURSE KRATOS_MPI_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

SET_SOURCE_FILES_PROPERTIES(${KRATOS_MPI_SEARCH_APPLICATION_SOURCES} ${KRATOS_MPI_APPLICATION_TESTING_SOURCES} PROPERTIES COMPILE_FLAGS -fno-strict-aliasing)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMPISearchApplication SHARED ${KRATOS_MPI_SEARCH_APPLICATION_SOURCES})
target_link_libraries(KratosMPISearchApplication KratosDEMApplication KratosCore ${MPI_LIBRARIES})
set_target_properties(KratosMPISearchApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMPISearchApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMPISearchApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosMPISearchApplication DESTINATION libs )


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MPISearchApplication.py" DESTINATION KratosMultiphysics )
