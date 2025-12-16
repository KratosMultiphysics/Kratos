set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosHDF5Application ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_HDF5_APPLICATION_CORE
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/hdf5_application_variables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_file.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_nd_data_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_points_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_connectivities_data.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_data_value_container_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_container_component_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/factor_elements_and_conditions_utility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/hdf5_data_set_partition_utility.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_properties_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_nodal_solution_step_variables_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_model_part_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_vertex_container_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_tensor_adaptor_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/hdf5_xdmf_connectivities_writer_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/vertex.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/vertex_utilities.cpp
)

## HDF5Application python interface sources
file(GLOB_RECURSE KRATOS_HDF5_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## HDF5Application testing
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_HDF5_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
  kratos_add_gtests(TARGET KratosHDF5Core SOURCES "${KRATOS_HDF5_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_definitions(-DH5_USE_110_API)

set(HDF5_USE_STATIC_LIBRARIES OFF)
if(USE_MPI)
  message(STATUS "Searching for parallel HDF5 libraries.")
  set(HDF5_PREFER_PARALLEL ON)
  set(KRATOS_HDF5_APPLICATION_CORE ${KRATOS_HDF5_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/hdf5_partitioned_model_part_io.cpp)
else(USE_MPI)
  message(STATUS "Searching for serial HDF5 libraries.")
endif(USE_MPI)

find_package(HDF5 REQUIRED COMPONENTS C)

if(HDF5_VERSION VERSION_LESS "1.8" )
  message(FATAL_ERROR "HDF5 1.8 is required.")
endif()

set(KRATOS_HDF5_APPLICATION_EXTRA_LIBS "")

if(USE_MPI)
  if(HDF5_IS_PARALLEL)
    message(STATUS "Using parallel HDF5 libraries.")
    if (   (HDF5_INCLUDE_DIRS MATCHES ".*openmpi.*")
        OR (HDF5_LIBRARY_DIRS MATCHES ".*openmpi.*"))
      add_definitions(-DOMPI_SKIP_MPICXX)
    endif()
    set(KRATOS_HDF5_APPLICATION_EXTRA_LIBS ${KRATOS_HDF5_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES})
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
    message(FATAL_ERROR "This is a parallel compilation. Parallel HDF5 libraries were not found.")
  endif()
else(USE_MPI)
  if(HDF5_IS_PARALLEL)
    message(FATAL_ERROR "This is a non-parallel compilation. Parallel HDF5 libraries were found.")
  else()
    message(STATUS "Using serial HDF5 libraries.")
  endif()
endif(USE_MPI)

include_directories(${HDF5_INCLUDE_DIRS})
add_definitions(${HDF5_DEFINITIONS})

add_library(KratosHDF5Core SHARED ${KRATOS_HDF5_APPLICATION_CORE})
if(USE_MPI)
  target_link_libraries(KratosHDF5Core PUBLIC KratosCore KratosMPICore ${HDF5_C_LIBRARIES} ${KRATOS_HDF5_APPLICATION_EXTRA_LIBS})
else(USE_MPI)
  target_link_libraries(KratosHDF5Core PUBLIC KratosCore ${HDF5_C_LIBRARIES} ${KRATOS_HDF5_APPLICATION_EXTRA_LIBS})
endif(USE_MPI)
set_target_properties(KratosHDF5Core PROPERTIES COMPILE_DEFINITIONS "HDF5_APPLICATION=EXPORT,API")

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosHDF5Application MODULE THIN_LTO ${KRATOS_HDF5_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosHDF5Application PRIVATE KratosHDF5Core)
set_target_properties(KratosHDF5Application PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosHDF5Core PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosHDF5Application PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosHDF5Application PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosHDF5Application PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/HDF5Application.py KratosMultiphysics/HDF5Application/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Fix Windows paths
file(TO_CMAKE_PATH "${HDF5_ROOT}/lib" HDF5_ROOT_DEP_LIB)
file(TO_CMAKE_PATH "${HDF5_ROOT}/bin" HDF5_ROOT_DEP_BIN)

# Install targets
install(TARGETS KratosHDF5Core DESTINATION libs)
install(TARGETS KratosHDF5Application DESTINATION libs)

# Install dependent libraries
# Note: Please be carefull to enable this option if MED or HDF5 are system wide dependencies as it will
# the full /lib system dir into kratos.
if(${INSTALL_TPL_LIBRARIES} MATCHES ON)
  if(WIN32)
      install(DIRECTORY "${HDF5_ROOT_DEP_LIB}/" DESTINATION libs FILES_MATCHING PATTERN "*.dll")
      install(DIRECTORY "${HDF5_ROOT_DEP_BIN}/" DESTINATION libs FILES_MATCHING PATTERN "*.dll")
  endif()
  if(APPLE)
      install(DIRECTORY "${HDF5_ROOT_DEP_LIB}/" DESTINATION libs FILES_MATCHING PATTERN "*.dylib")
      install(DIRECTORY "${HDF5_ROOT_DEP_BIN}/" DESTINATION libs FILES_MATCHING PATTERN "*.dylib")
  endif()
  if(UNIX AND NOT APPLE)
      file(GLOB_RECURSE LIB_FILES "${HDF5_ROOT_DEP_LIB}/*.so*")
      file(GLOB_RECURSE LIB_FILES APPEND "${HDF5_ROOT_DEP_BIN}/*.so*")
      install(FILES ${LIB_FILES} DESTINATION libs)
  endif()
endif()

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosHDF5Core" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosHDF5Application" PARENT_SCOPE)