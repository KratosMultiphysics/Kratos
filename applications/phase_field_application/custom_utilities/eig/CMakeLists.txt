CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#find_package(Boost REQUIRED)
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}) 
#ADD_DEFINITIONS(-g -O2 -fsigned-char -freg-struct-return -Wall -W -Wshadow -Wstrict-prototypes -Wpointer-arith -Wcast-qual -Winline -fpermissive)
ADD_DEFINITIONS(-g -O2)

add_library(eig STATIC eig3.cpp)

set(file_list
test
)

foreach(str ${file_list})
    ADD_EXECUTABLE(${str} ${str}.cpp)
    TARGET_LINK_LIBRARIES(${str} eig)
endforeach()

