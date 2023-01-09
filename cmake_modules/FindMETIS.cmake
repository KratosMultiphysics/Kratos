# Find the METIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# METIS_INCLUDE_DIR - where to find autopack.h
# METIS_LIBRARIES   - List of fully qualified libraries to link against.
# METIS_FOUND       - Do not attempt to use if "no" or undefined.

# Manual setting if the METIS_ROOT_DIR is set
IF (METIS_ROOT_DIR)
    SET(METIS_INCLUDE_DIR "${METIS_ROOT_DIR}/include")
    SET(METIS_LIBRARY_DIR "${METIS_ROOT_DIR}/lib")
ENDIF (METIS_ROOT_DIR)

FIND_PATH(METIS_INCLUDE_DIR metis.h
    /usr/local/include
    /usr/include
    "${METIS_INCLUDE_DIR}"
)

FIND_LIBRARY(METIS_LIBRARY metis
    /usr/local/lib
    /usr/lib
    "${METIS_LIBRARY_DIR}"
)

IF(METIS_INCLUDE_DIR AND METIS_LIBRARY)
    SET( METIS_LIBRARIES ${METIS_LIBRARY})
    SET( METIS_FOUND TRUE )
    MESSAGE(STATUS "METIS library found")
    MESSAGE(STATUS "METIS_INCLUDE_DIR: ${METIS_INCLUDE_DIR}")
    MESSAGE(STATUS "METIS_LIBRARIES: ${METIS_LIBRARIES}")
ELSE(METIS_INCLUDE_DIR AND METIS_LIBRARY)
    MESSAGE(FATAL_ERROR "Finding METIS failed, please try to set the var METIS_ROOT_DIR, or METIS_INCLUDE_DIR and METIS_LIBRARY_DIR")
ENDIF(METIS_INCLUDE_DIR AND METIS_LIBRARY)
