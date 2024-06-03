set(CMAKE_INCLUDE_CURRENT_DIR ON)

# generate variables with the sources
file(GLOB_RECURSE QUESO_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/embedding/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/quadrature/integration_points_1d/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/containers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/includes/*.cpp

    # Add external libraries
    ${CMAKE_CURRENT_SOURCE_DIR}/../aabb_tree/*.cc
)

add_library(queso STATIC ${QUESO_SOURCES})

set(QUESO_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(QUESO_LIBRARIES queso)

message(STATUS "queso : ${QUESO_INCLUDE_DIR}" )
