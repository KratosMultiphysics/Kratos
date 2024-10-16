# G+Smo
find_package(gismo REQUIRED)
include_directories( ${GISMO_INCLUDE_DIRS})
link_directories   ( ${GISMO_LIBRARY_DIR} )

include( ${GISMO_USE_FILE} )
set(CMAKE_MODULE_PATH ${GISMO_MODULE_PATH})

# include(gsConfig)

message(STATUS "G+Smo shared library found in ${GISMO_INCLUDE_DIRS}")
