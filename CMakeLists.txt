cmake_minimum_required (VERSION 2.8)
project (KratosBaseProject)

#set compiler flags
if(${CMAKE_COMPILER_IS_GNUCXX}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -funroll-loops -ffast-math  ")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -funroll-loops -ffast-math  ")
  message("additional default options were set for gcc")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_COMPILER_IS_GNUCXX}) 
if(${CMAKE_C_COMPILER} MATCHES "icc.*$") 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -ansi -funroll-loops  -lpthread -wd654 -wd10010  ")
  set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fPIC -funroll-loops  -lpthread -wd654 -wd10010  ")
  message("additional default options were set for intel compiler")
  message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")
  message("CMAKE_C_FLAGS = ${CMAKE_C_FLAGS}")
endif(${CMAKE_C_COMPILER} MATCHES "icc.*$") 


#tell the lionker to give an error if undefined functions are found
set(CMAKE_SHARED_LINKER_FLAGS "-Wl,--no-undefined")

#set kratos specific module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

######################################################################################
######################################################################################
######################################################################################
##find all dependencies

##*****************************
#find and include OpenMP if possible
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif (OPENMP_FOUND) 

##*****************************
#finding and including BOOST library
set(Boost_USE_STATIC_LIBS   OFF)
set(Boost_USE_MULTITHREADED ON)
find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

##*****************************
#finding and including PYTHON
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

##*****************************
#finding blas and lapack if needed
set(BLAS_INCLUDE_NEEDED OFF)
if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )
    set(BLAS_INCLUDE_NEEDED ON)
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON )

if(${BLAS_INCLUDE_NEEDED} MATCHES ON )
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
	enable_language( Fortran )
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

    find_package(BLAS REQUIRED)
    include_directories(${BLAS_INCLUDE_DIRS})
    link_directories(${BLAS_LIBRARIES})

    find_package(LAPACK REQUIRED)
    include_directories(${LAPACK_INCLUDE_DIRS})
    link_directories(${LAPACK_LIBRARIES})
endif(${BLAS_INCLUDE_NEEDED} MATCHES ON )

##*****************************
#finding Mpi
set(MPI_NEEDED OFF)
if(${METIS_APPLICATION} MATCHES ON )
    set(MPI_NEEDED ON)
endif(${METIS_APPLICATION} MATCHES ON )
if(${TRILINOS_APPLICATION} MATCHES ON )
    set(MPI_NEEDED ON)
endif(${TRILINOS_APPLICATION} MATCHES ON )

if(${MPI_NEEDED} MATCHES ON )
    find_package(MPI REQUIRED)
    add_definitions( -DKRATOS_USING_MPI )
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_C_COMPILE_FLAGS ${CMAKE_C_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
    set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

    #note that this will be included in the compilation of the kratos core, so that it is not needed to include it somewhere else
    include_directories(${MPI_INCLUDE_PATH})
endif(${MPI_NEEDED} MATCHES ON )


##echo user options
message( " ")
message("CMAKE_BUILD_TYPE ........ ${CMAKE_BUILD_TYPE}")
message( " ")
######################################################################################
######################################################################################
######################################################################################
#include internal dependencies
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/omptl )
# include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/gidpost )

#defines needed
add_definitions( -DKRATOS_PYTHON )
add_definitions( -DNDEBUG )
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_definitions( -fPIC )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#include subdirectories
add_subdirectory(external_libraries/gidpost)
add_subdirectory(kratos)
add_subdirectory(applications)

#set installation directory
set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR} )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("inc_dirs = ${inc_dirs}")