# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR - where to find autopack.h
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against.
# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
  /usr/local/include
  /usr/include
  "${PARMETIS_ROOT_DIR}"
)

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
  /usr/local/lib
  /usr/lib
  "${PARMETIS_ROOT_DIR}"
)

FIND_LIBRARY(METIS_LIBRARY metis
  /usr/local/lib
  /usr/lib
  "${PARMETIS_ROOT_DIR}"
)

IF(PARMETIS_INCLUDE_DIR)
  IF(PARMETIS_LIBRARY)
    SET( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
    SET( PARMETIS_FOUND "YES" )
    message("PARMETIS found")
  ELSE(PARMETIS_LIBRARY)
    message("finding PARMETIS failed, please try to set the var PARMETIS_ROOT_DIR")
  ENDIF(PARMETIS_LIBRARY)
ELSE(PARMETIS_INCLUDE_DIR)
  message("finding PARMETIS failed, please try to set the var PARMETIS_ROOT_DIR")
ENDIF(PARMETIS_INCLUDE_DIR)
