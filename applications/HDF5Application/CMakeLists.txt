set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosHDF5Application ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_HDF5_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_io.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_file.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/hdf5_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/tests/test_hdf5_file.cpp  
)

set(HDF5_USE_STATIC_LIBRARIES OFF)
if(MPI_NEEDED)
  message(STATUS "Searching for parallel HDF5 libraries.")
  set(HDF5_PREFER_PARALLEL ON)
  set(KRATOS_HDF5_APPLICATION_SOURCES ${KRATOS_HDF5_APPLICATION_SOURCES} ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_file_mpi.cpp)
else()
  message(STATUS "Searching for serial HDF5 libraries.")
endif()

find_package(HDF5 REQUIRED COMPONENTS C)

string(REGEX MATCH "1\\.[0-9]+" HDF5_MAJOR_MINOR_VERSION ${HDF5_VERSION})
if(NOT HDF5_MAJOR_MINOR_VERSION EQUAL "1.8")
  message(FATAL_ERROR "HDF5 1.8 is required.")
endif()

if(MPI_NEEDED)
  if(HDF5_IS_PARALLEL)
    message(STATUS "Using parallel HDF5 libraries.")
    if (   (HDF5_INCLUDE_DIRS MATCHES ".*openmpi.*")
        OR (HDF5_LIBRARY_DIRS MATCHES ".*openmpi.*"))
      add_definitions(-DOMPI_SKIP_MPICXX)
    endif()
    ##################################################################
    # If cmake configures another application using HDF5 in serial,
    # the find_library() used in FindHDF5.cmake may find the serial
    # version of libhdf5.so even though the directory passed with
    # HINTS is for the parallel version. If the other application
    # cannot be switched off, a workaround is to call find_library
    # with NO_CMAKE_PATH:
    #    find_library(HDF5_LIBRARY
    #                 NAMES "hdf5"
    #                 HINTS ${HDF5_C_LIBRARY_DIRS}
    #                 ENV HDF5_ROOT
    #                 PATH_SUFFIXES lib Lib 
    #                 NO_CMAKE_PATH)
    ##################################################################
  else()
    message(FATAL_ERROR "Parallel HDF5 libraries were not found.")
  endif()
else()
  if(HDF5_IS_PARALLEL)
    message(FATAL_ERROR "Parallel HDF5 libraries were found.")
  else()
    message(STATUS "Using serial HDF5 libraries.")
  endif()
endif()

include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})

# define library Kratos which defines the basic python interface
add_library(KratosHDF5Application SHARED ${KRATOS_HDF5_APPLICATION_SOURCES})
target_link_libraries(KratosHDF5Application KratosCore ${HDF5_C_LIBRARIES})
set_target_properties(KratosHDF5Application PROPERTIES PREFIX "")
install(TARGETS KratosHDF5Application DESTINATION libs)

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
