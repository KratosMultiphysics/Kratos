
# include(FindPackageHandleStandardArgs)
# find_package_handle_standard_args






# Find the MED includes and libraries
#
# MED_INCLUDE_DIR - where to find autopack.h
# MED_LIBRARIES   - List of fully qualified libraries to link against.
# MED_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(MED_INCLUDE_DIR "med.h"
  PATHS
  "${MED_ROOT}/include"
  /usr/local/include
  /usr/include
)

FIND_LIBRARY(MED_LIBRARY medC
  PATHS
  "${MED_ROOT}/lib"
  /usr/local/lib/*
  /usr/lib/*
)

IF(MED_INCLUDE_DIR AND MED_LIBRARY)
  SET( MED_LIBRARIES ${MED_LIBRARY})
  SET( MED_FOUND "YES" )
  message(STATUS "MED found!")
  message(STATUS "MED includes: ${MED_INCLUDE_DIR}")
  message(STATUS "MED libraries: ${MED_LIBRARY}")
else()
  message(WARNING "finding MED failed, please try to set the var MED_ROOT")
ENDIF()
