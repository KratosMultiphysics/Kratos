# Not including current dir (we want to force users to specify the mpi/ in the include path)
set(CMAKE_INCLUDE_CURRENT_DIR OFF)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## Kratos MPI extension sources
set( KRATOS_MPI_SOURCES
   mpi_environment.cpp
   sources/mpi_data_communicator.cpp
)

## Kratos MPI extension Python interface
file(GLOB_RECURSE KRATOS_MPI_PYTHON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp)

## Kratos MPI test sources.
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MPI_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/sources/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Define library for the C++ layer
add_library(KratosMPICore SHARED ${KRATOS_MPI_SOURCES} ${KRATOS_MPI_TEST_SOURCES})
target_link_libraries(KratosMPICore PUBLIC KratosCore ${MPI_LIBRARIES} )
#set_target_properties(KratosMPICore PROPERTIES COMPILE_DEFINITIONS "KRATOS_MPI_CORE=IMPORT,API")

## Define library for the Python interface layer
pybind11_add_module(KratosMPI MODULE THIN_LTO ${KRATOS_MPI_PYTHON_SOURCES})
# add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(KratosMPI PRIVATE KratosMPICore KratosCore )
# set_target_properties(KratosMPI PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosMPI)
endif(USE_COTIRE MATCHES ON)

# Install library
install(TARGETS KratosMPICore DESTINATION libs )
install(TARGETS KratosMPI DESTINATION libs)

# Install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/mpi_module_init.py" DESTINATION "KratosMultiphysics/mpi" RENAME "__init__.py" )
