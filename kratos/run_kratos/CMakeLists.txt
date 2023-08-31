set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( EMBEDDED_PYTHON_SOURCES
 	${CMAKE_CURRENT_SOURCE_DIR}/krun_main.cpp
     )


###############################################################
message("python libraries to be installed are: " ${PYTHON_LIBRARIES} )
install(FILES ${PYTHON_LIBRARIES} DESTINATION .)
