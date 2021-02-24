set(CMAKE_INCLUDE_CURRENT_DIR ON)

# generate variables with the sources
set(CLIPPER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/clipper.cpp)

add_library(clipper STATIC ${CLIPPER_SOURCES})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries(clipper m)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

set(CLIPPER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CLIPPER_LIBRARIES clipper)

message(STATUS "clipper : ${CLIPPER_INCLUDE_DIR}" )
# install(TARGETS triangle DESTINATION libs )
