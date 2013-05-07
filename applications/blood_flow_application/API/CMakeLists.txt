set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(${PYTHON_INCLUDE_DIRS})
link_directories(${PYTHON_LIBRARIES})

############ module to include sys libs #############
include(InstallRequiredSystemLibraries)

## generate variables with the sources
set( CVP_LIB_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/cvp.cpp
     ) 

#if(WIN32)
#	include(GenerateExportHeader)
#	add_compiler_export_flags()
#endif()

###############################################################
add_library(cvp SHARED ${CVP_LIB_SOURCES} )
# target_link_libraries(cvp C:/Python27/libs/python27.lib)

target_link_libraries(cvp ${EXTRA_LINK_LIBS} ${PYTHON_LIBRARIES})

if(WIN32)
##generate_export_header(cvp)
##	GENERATE_EXPORT_HEADER( cvp 
##		BASE_NAME cvp
##		EXPORT_MACRO_NAME cvp_EXPORT
##		EXPORT_FILE_NAME cvp_Export.h
##		STATIC_DEFINE cvp_BUILT_AS_STATIC
##		)
	install(FILES cvp.h DESTINATION libs)
	install(TARGETS cvp DESTINATION libs )
else()
	install(TARGETS cvp LIBRARY DESTINATION libs )
endif()	

add_executable(testcvp ${CMAKE_CURRENT_SOURCE_DIR}/testcvp.cpp )
target_link_libraries(testcvp cvp )
install(TARGETS testcvp DESTINATION libs )

message("python libraries to be installed are: " ${PYTHON_LIBRARIES} )
install(FILES ${PYTHON_LIBRARIES} DESTINATION .)

