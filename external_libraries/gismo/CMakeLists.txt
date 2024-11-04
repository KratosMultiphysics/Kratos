# G+Smo
find_package(gismo REQUIRED)

include( ${GISMO_USE_FILE} )

#set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})

set(GISMO_LIBRARY_DIR ${GISMO_LIBRARY_DIR} PARENT_SCOPE)
set(GISMO_INCLUDE_DIRS ${GISMO_INCLUDE_DIRS} PARENT_SCOPE)
link_directories ( ${GISMO_LIBRARY_DIR} )

 

# include(gsConfig)

message(STATUS "G+Smo shared library found in ${GISMO_LIBRARY_DIR}")

message(STATUS "G+Smo include directories ${GISMO_INCLUDE_DIRS}")
