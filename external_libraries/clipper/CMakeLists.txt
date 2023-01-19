set(CMAKE_INCLUDE_CURRENT_DIR ON)

# generate variables with the sources
set(CLIPPER_SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/src/clipper.engine.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/clipper.offset.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/clipper.rectclip.cpp
)

add_library(clipper STATIC ${CLIPPER_SOURCES})

set(CLIPPER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CLIPPER_LIBRARIES clipper)

message(STATUS "clipper : ${CLIPPER_INCLUDE_DIR}" )
