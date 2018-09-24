
include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## Kratos MPI extension sources
set( KRATOS_MPI_SOURCES
   dummy.cpp
)

## Kratos MPI extension Python interface
set( KRATOS_MPI_PYTHON_SOURCES
    python/kratos_mpi_python.cpp
)

## Define library for the C++ layer
add_library(KratosMPICore SHARED ${KRATOS_MPI_SOURCES} )
target_link_libraries(KratosMPICore PUBLIC KratosCore ${MPI_LIBRARIES} )
#set_target_properties(KratosMPICore PROPERTIES COMPILE_DEFINITIONS "KRATOS_MPI_CORE=IMPORT,API")

## Define library for the Python interface layer
pybind11_add_module(KratosMPI MODULE THIN_LTO ${KRATOS_MPI_PYTHON_SOURCES})
# add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(KratosMPI PRIVATE KratosMPICore )
# set_target_properties(KratosMPI PROPERTIES PREFIX "")

# Install library
install(TARGETS KratosMPICore DESTINATION libs )
install(TARGETS KratosMPI DESTINATION libs)

# Install python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/KratosMultiphysicsMPI.py" DESTINATION KratosMultiphysics )
