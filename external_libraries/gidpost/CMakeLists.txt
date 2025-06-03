cmake_minimum_required (VERSION 2.9)

project (gidpost)

## # Add your project files and targets here
## add_executable(MyExecutable main.cpp)
## 
## # Link against the found libraries
## target_link_libraries(MyExecutable ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES})

# The gidpost version number,
# GP_VERSION_MAJOR, GP_VERSION_MINOR and GP_VERSION are defined in gidpost.h
# set (GP_VERSION_MAJOR 2)
# set (GP_VERSION_MINOR 10)
# set( GP_VERSION ${GP_VERSION_MAJOR}.${GP_VERSION_MINOR} )

# Find zlib
find_package(ZLIB REQUIRED)
if(ZLIB_FOUND)
    message(STATUS "Found ZLIB: ${ZLIB_INCLUDE_DIRS}")
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

# Find hdf5
find_package(HDF5 REQUIRED COMPONENTS C HL)
if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_INCLUDE_DIRS}")
    include_directories(${HDF5_INCLUDE_DIRS})
    # HDF5 variable used in configure_file( ... gidpost_config.h.in)
    # set (HDF5 1)
    option (ENABLE_HDF5 "Use HDF5" ON)
else(HDF5_FOUND)
    option (ENABLE_HDF5 "Use HDF5" OFF)
endif(HDF5_FOUND)

option(ENABLE_SHARED_LIBS "Build gidpost as a shared libraries." OFF)

if (ENABLE_SHARED_LIBS)
  add_definitions (-DGIDPOST_SHARED)
endif (ENABLE_SHARED_LIBS)
option ( ENABLE_FORTRAN_EXAMPLES "Build fortran examples, it is needed a fortran compiler" OFF )

if (ENABLE_FORTRAN_EXAMPLES)
  enable_language ( Fortran )
  # add_definitions(-DgFortran)
  # message(STATUS "*** defining -DgFortran")
  add_definitions (-DgFortran)
endif (ENABLE_FORTRAN_EXAMPLES)

# used in examples/MultiplePartition
option(ENABLE_PARALLEL_EXAMPLE "Build gidpost as a shared libraries." OFF)
if( ENABLE_PARALLEL_EXAMPLE )
  find_package( Boost COMPONENTS thread REQUIRED)
endif( ENABLE_PARALLEL_EXAMPLE )

# using g++ 11 in older linuxes
# to enable debug information
if (CMAKE_COMPILER_IS_GNUCXX)
  # get kernel version
  execute_process( COMMAND uname -r OUTPUT_VARIABLE UNAME_RESULT OUTPUT_STRIP_TRAILING_WHITESPACE)
  message( -- " Kernel version: " ${UNAME_RESULT})
  string( REGEX MATCH "[0-9]+.[0-9]+" LINUX_KERNEL_VERSION ${UNAME_RESULT})
  # if ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0)
  if ( LINUX_KERNEL_VERSION VERSION_LESS 5.0)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -gdwarf-4")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -gdwarf-4")
    message(STATUS "NEW FLAGS = ${CMAKE_CXX_FLAGS}")
  # endif ( CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 11.0)
  endif ( LINUX_KERNEL_VERSION VERSION_LESS 5.0)
endif (CMAKE_COMPILER_IS_GNUCXX)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/gidpost_config.h.in"
  "${PROJECT_BINARY_DIR}/gidpost_config.h"
)

# so that gidpost library uses the above gidpost_config.h file
add_definitions (-DHAVE_GIDPOST_CONFIG_H)

install( FILES "${PROJECT_BINARY_DIR}/gidpost_config.h" DESTINATION "include" )

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

include_directories ("${PROJECT_SOURCE_DIR}")
include_directories ("${PROJECT_SOURCE_DIR}/source")

add_subdirectory (source)

# Make sure the linker can find the gidpost library once it is
# built.
link_directories (${PROJECT_BINARY_DIR}/source  ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES})

add_subdirectory( examples )
add_subdirectory( examples/MultiplePartition )
