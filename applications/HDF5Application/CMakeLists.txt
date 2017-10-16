set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosHDF5Application ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_HDF5_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_io.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/hdf5_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
)

set( HDF5_IO ON )
if( ${HDF5_IO} MATCHES ON )
  set( HDF5_USE_STATIC_LIBRARIES OFF )
  if( ${MPI_NEEDED} MATCHES ON )
    set( HDF5_PREFER_PARALLEL ON )
  endif( ${MPI_NEEDED} MATCHES ON )
#  find_package(HDF5 CXX)
  find_package(HDF5 COMPONENTS C REQUIRED)

  include_directories( ${HDF5_INCLUDE_DIRS} )
endif( ${HDF5_IO} MATCHES ON )

# define library Kratos which defines the basic python interface
add_library(KratosHDF5Application SHARED ${KRATOS_HDF5_APPLICATION_SOURCES})
target_link_libraries(KratosHDF5Application KratosCore ${HDF5_C_LIBRARIES} )
set_target_properties(KratosHDF5Application PROPERTIES PREFIX "")
install(TARGETS KratosHDF5Application DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosHDF5Application PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosHDF5Application PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/HDF5Application.py" DESTINATION KratosMultiphysics )
