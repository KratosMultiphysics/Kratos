cmake_minimum_required (VERSION 2.8)
project (KratosBaseProject)

##echo user options
message( " ")
message("CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}")
message( " ")

######################################################################################
##find all dependencies

#find and include OpenMP if possible
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OPENMP_FOUND) 

#finding and including BOOST library
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

#finding and including PYTHON
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})


enable_language( Fortran )
#finding blas
find_package(BLAS)

#finding lapack
find_package(LAPACK)

#finding Mpi
# find_package(Trilinos) 



######################################################################################
#include internal dependencies
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/omptl )
# include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/gidpost )

#defines needed
add_definitions( -DKRATOS_PYTHON )
add_definitions( -DNDEBUG )

#include subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)
add_subdirectory(applications)

#set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("inc_dirs = ${inc_dirs}")