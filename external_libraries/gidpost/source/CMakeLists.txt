set (gidpost_source 
  gidpost.c gidpostInt.c  gidpostfor.c gidpostforAPI.c
  gidpostHash.c hashtab.c recycle.c lookupa.c)

add_definitions(-DNDEBUG)
message("WARNING: YOUR GIDPOST IS BEING COMPILED  WITHOUT DEBUG (ALL ASSERTS WILL BE AVOIDED)")

find_package(HDF5)
if (HDF5_FOUND)
  option (HDF5 "Use HDF5" OFF)
endif(HDF5_FOUND)

if (HDF5)
  set (hdf5_source hdf5c.c gidpostHDF5.c)
endif (HDF5)

if (BUILD_SHARED_LIBS)
  set (LIB_TYPE SHARED)
else (BUILD_SHARED_LIBS)
  set (LIB_TYPE STATIC)
endif (BUILD_SHARED_LIBS)

add_library (gidpost ${LIB_TYPE} ${gidpost_source} ${hdf5_source})

if (BUILD_SHARED_LIBS)
  add_definitions (-DGIDPOST_EXPORTS)
endif (BUILD_SHARED_LIBS)

if (MSVC)
  add_definitions (-D_CRT_SECURE_NO_DEPRECATE)
endif (MSVC)

if (HDF5)
  include_directories(${HDF5_INCLUDE_DIRS})
  target_link_libraries(gidpost ${HDF5_LIBRARIES})
endif (HDF5)

find_package(ZLIB REQUIRED)
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

install (FILES gidpost.h ${gidpost_BINARY_DIR}/gidpost_config.h
  DESTINATION include)
