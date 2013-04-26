set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( CVP_LIB_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/cvp.cpp
     ) 


###############################################################
add_library(cvp ${CVP_LIB_SOURCES} )
target_link_libraries(cvp ${EXTRA_LINK_LIBS} ${PYTHON_LIBRARIES})
#install(TARGETS cvp DESTINATION libs )

add_executable(testcvp ${CMAKE_CURRENT_SOURCE_DIR}/testcvp.cpp )
target_link_libraries(testcvp cvp )
install(TARGETS testcvp DESTINATION libs )

message("python libraries to be installed are: " ${PYTHON_LIBRARIES} )
install(FILES ${PYTHON_LIBRARIES} DESTINATION .)

