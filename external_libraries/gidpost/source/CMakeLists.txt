#include (ExternalProject)
set (gidpost_source 
  gidpost.c gidpost_cluster.c gidpostFILES.c gidpostInt.c  gidpostfor.c gidpostforAPI.c
  gidpostHash.c hashtab.c recycle.c lookupa.c)

  
# Remove Unity ofr gidpost
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/gidpost.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/gidpostInt.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/gidpostHash.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/hashtab.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/recycle.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
# set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/source/lookupa.c PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)

# find  package already done in  ../CmakeLists.txt
# find_package(HDF5)
# if (HDF5_FOUND)
#   option (ADD_HDF5_SUPPORT "Use HDF5" ON)
# endif(HDF5_FOUND)

if (ENABLE_HDF5)
  set (hdf5_source hdf5c.c gidpostHDF5.c)
endif (ENABLE_HDF5)

if (ENABLE_SHARED_LIBS)
  set (LIB_TYPE SHARED)
else (ENABLE_SHARED_LIBS)
  set (LIB_TYPE STATIC)
endif (ENABLE_SHARED_LIBS)

add_library (gidpost ${LIB_TYPE} ${gidpost_source} ${hdf5_source})

if (ENABLE_SHARED_LIBS)
  add_definitions (-DGIDPOST_EXPORTS)
endif (ENABLE_SHARED_LIBS)

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_DEPRECATE)
endif (MSVC)

# if()
#   set(CMAKE_CXX_FLAGS_DEBUG_INIT "-Wall")
# endif()
if ( ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU") OR ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang"))
  target_compile_options( gidpost PRIVATE "$<$<CONFIG:Debug>:-Wall>")
  target_compile_options( gidpost PRIVATE "$<$<CONFIG:Debug>:-Wextra>")
  target_compile_options( gidpost PRIVATE "$<$<CONFIG:Debug>:-Wconversion>")
  # add_compile_options( "$<$<CONFIG:Debug>:-Wall>")
  # add_compile_options( "$<$<CONFIG:Debug>:-Wextra>")
endif()

if (ENABLE_HDF5)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(gidpost ${HDF5_LIBRARIES})
endif (ENABLE_HDF5)

# find  package already done in  ../CmakeLists.txt
# find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
  target_link_libraries(gidpost ${ZLIB_LIBRARIES})
else ( ZLIB_FOUND )
  message ( SEND_ERROR "ZLIB_LIBRARY NOTFOUND, set ZLIB_INCLUDE_DIRS and ZLIB_LIBRARIES variables")
endif(ZLIB_FOUND)

find_package( Threads )
if ( Threads_FOUND )
  target_link_libraries(gidpost ${CMAKE_THREAD_LIBS_INIT})
endif (Threads_FOUND)

if (UNIX)
  target_link_libraries(gidpost m)
endif (UNIX)

# install rules
#
install (TARGETS gidpost
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

install( FILES gidpost.h gidpost_types.h gidpost_functions.h gidpost_cluster_functions.h gidpost_functions_deprecated.h DESTINATION include )
