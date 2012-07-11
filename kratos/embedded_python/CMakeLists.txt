set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

## generate variables with the sources
set( EMBEDDED_PYTHON_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/krun_main.cpp
     )


###############################################################
add_executable(krun ${EMBEDDED_PYTHON_SOURCES} )
target_link_libraries(krun  ${PYTHON_LIBRARIES})
install(TARGETS krun DESTINATION . )
install(FILES ${PYTHON_LIBRARIES} DESTINATION .)