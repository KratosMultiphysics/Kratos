# cmake_minimum_required (VERSION 2.8)
# project (MultiplePartition)

# find_path( GIDPOST_INCLUDE_DIR
#   NAMES gidpost.h
#   PATHS ${GIDPOST_PREFIX_DIR}/include )
# 
# if( GIDPOST_INCLUDE_DIR )
#   message( "GIDPOST_INCLUDE_DIR = ${GIDPOST_INCLUDE_DIR}" )
#   include_directories( "${GIDPOST_INCLUDE_DIR}" )
# else( GIDPOST_INCLUDE_DIR )
#   message( FATAL_ERROR "Unable to find gidpost.h, set GIDPOST_PREFIX_DIR" )
# endif( GIDPOST_INCLUDE_DIR )
# 
# find_library( GIDPOST_LIBRARY
#   NAMES gidpost
#   PATHS ${GIDPOST_PREFIX_DIR}/lib )
# 
# if( GIDPOST_LIBRARY )
#   message( "GIDPOST_LIBRARY = ${GIDPOST_LIBRARY}" )
#   list( APPEND GIDPOST_LIBRARIES ${GIDPOST_LIBRARY} )
# else( GIDPOST_LIBRARY )
#   message( FATAL_ERROR "Unable to find library gidpost, set GIDPOST_PREFIX_DIR" )
# endif( GIDPOST_LIBRARY )

include_directories ("${PROJECT_SOURCE_DIR}/examples")

add_executable( SimulateMultiplePartition SimulateMultiplePartition.cxx )

# target_link_libraries( SimulateMultiplePartition ${GIDPOST_LIBRARIES} )
target_link_libraries( SimulateMultiplePartition gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

# Find Boost
# find_package(Boost REQUIRED)

if( ENABLE_PARALLEL_EXAMPLE)
  if( Boost_FOUND )
    include_directories( ${Boost_INCLUDE_DIRS} )
    add_definitions( -DHAVE_BOOST_THREAD )
    target_link_libraries( SimulateMultiplePartition ${Boost_LIBRARIES} )
  else( Boost_FOUND )
    message( WARNING "Boost thread component not found, parallel feature will be disabled" )
  endif( Boost_FOUND )
endif( ENABLE_PARALLEL_EXAMPLE)
