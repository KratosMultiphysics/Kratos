set(CMAKE_INCLUDE_CURRENT_DIR ON)
cmake_minimum_required(VERSION 3.9)

if(CSV_CXX_STANDARD)
	set(CMAKE_CXX_STANDARD ${CSV_CXX_STANDARD})
else()
	set(CMAKE_CXX_STANDARD 11)
endif(CSV_CXX_STANDARD)

message("Building CSV library using C++${CMAKE_CXX_STANDARD}")

# Defines CSV_HAS_CXX17 in compatibility.hpp
if (CMAKE_VERSION VERSION_LESS "3.12.0")
	add_definitions(-DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD})
else()
	add_compile_definitions(CMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD})
endif()

set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads QUIET REQUIRED)

set(CSV_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/include/)
set(CSV_SOURCE_DIR ${CSV_INCLUDE_DIR}/internal/)

include_directories(${CSV_INCLUDE_DIR})

## Main Library
add_subdirectory(${CSV_SOURCE_DIR})

###############################################################

add_library(csv_parser STATIC ${CSV_SOURCE_DIR})

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  target_link_libraries(csv_parser m)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

set(CSV_PARSER_LIBRARIES csv_parser)

message(STATUS "csv-parser : ${CSV_INCLUDE_DIR}")
