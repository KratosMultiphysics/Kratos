set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

message("WARNING: Deprecation: runkratos will be deleted by the end of 2022, please consider disabling it")

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( EMBEDDED_PYTHON_SOURCES
 	${CMAKE_CURRENT_SOURCE_DIR}/krun_main.cpp
     )


###############################################################
add_executable(runkratos ${EMBEDDED_PYTHON_SOURCES} )
target_link_libraries(runkratos ${EXTRA_LINK_LIBS} ${PYTHON_LIBRARIES})

install(TARGETS runkratos DESTINATION . )

message("python libraries to be installed are: " ${PYTHON_LIBRARIES} )
install(FILES ${PYTHON_LIBRARIES} DESTINATION .)
