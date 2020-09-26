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

include_directories(${CSV_INCLUDE_DIR})
