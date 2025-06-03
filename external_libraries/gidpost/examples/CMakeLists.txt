include_directories ("${PROJECT_SOURCE_DIR}/examples")

# deprecated API example:
add_executable ( testc testpost.c )
target_link_libraries ( testc gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

add_executable ( ex1_sphere_fd ex1_sphere_fd.c )
target_link_libraries ( ex1_sphere_fd gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

add_executable ( IGA_test_1 IGA_test_1.c )
target_link_libraries ( IGA_test_1 gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

## BugMesh* are sources to do some debugging
## add_executable ( BugMeshGroup BugMeshGroup.c )
## target_link_libraries ( BugMeshGroup gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
## 
## add_executable ( BugMeshGroup02 BugMeshGroup02.c )
## target_link_libraries ( BugMeshGroup02 gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
## 
## add_executable ( BugMeshGroup03 BugMeshGroup03.c )
## target_link_libraries ( BugMeshGroup03 gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )

if (CMAKE_COMPILER_IS_GNUCXX)
  # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -D_BSD_SOURCE")
  # _POSIX_C_SOURCE needed to get strdup() declaration as it is not std c99
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -D_POSIX_C_SOURCE=200809L")
endif (CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_Fortran_COMPILER)
  add_executable ( testf testpostfor.f )
  target_link_libraries ( testf gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
  
  add_executable ( test_fd_f90 testpost_fd.f90 ../fortran_module/gidpost_interface.f90)
  target_link_libraries ( test_fd_f90 gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
  
  if (UNIX)
    target_link_libraries ( testf gfortran )
    target_link_libraries ( test_fd_f90 gfortran ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
  endif (UNIX)
endif (CMAKE_Fortran_COMPILER)

add_executable( TestMultiplePartition TestMultiplePartition.cxx )
target_link_libraries( TestMultiplePartition gidpost ${ZLIB_LIBRARIES} ${HDF5_LIBRARIES} ${HDF5_HL_LIBRARIES} )
