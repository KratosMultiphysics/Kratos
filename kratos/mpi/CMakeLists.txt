# Not including current dir (we want to force users to specify the mpi/ in the include path)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

message("**** configuring Kratos to be compiled with MPI-Support ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Kratos MPI extension sources
set( KRATOS_MPI_SOURCES
   sources/mpi_data_communicator.cpp
   sources/mpi_manager.cpp
   utilities/data_communicator_factory.cpp
   utilities/gather_modelpart_utility.cpp
   utilities/mpi_normal_calculation_utilities.cpp;
   utilities/parallel_fill_communicator.cpp
)

## Kratos MPI extension Python interface
file(GLOB_RECURSE KRATOS_MPI_PYTHON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

## Kratos MPI test sources.
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MPI_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/sources/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################

## Define library for the C++ layer
add_library(KratosMPICore SHARED ${KRATOS_MPI_SOURCES} ${KRATOS_MPI_TEST_SOURCES})
target_link_libraries(KratosMPICore PUBLIC KratosCore ${MPI_LIBRARIES} )
#set_target_properties(KratosMPICore PROPERTIES COMPILE_DEFINITIONS "KRATOS_MPI_CORE=IMPORT,API")

## Define library for the Python interface layer
pybind11_add_module(KratosMPI MODULE THIN_LTO ${KRATOS_MPI_PYTHON_SOURCES})
# add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(KratosMPI PRIVATE KratosMPICore KratosCore )
# set_target_properties(KratosMPI PROPERTIES PREFIX "")


###############################################################################

if(USE_COTIRE MATCHES ON)
    cotire(KratosMPICore)
    cotire(KratosMPI)
endif(USE_COTIRE MATCHES ON)

###############################################################################

# Install library
install(TARGETS KratosMPICore DESTINATION libs )
install(TARGETS KratosMPI DESTINATION libs)

# Install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/mpi_module_init.py" DESTINATION "KratosMultiphysics/mpi" RENAME "__init__.py" )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

if(${INSTALL_TESTING_FILES} MATCHES ON )
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION kratos/mpi
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)
