set(CMAKE_INCLUDE_CURRENT_DIR ON)

## generate variables with the sources
set( BPRINTER_SOURCES     
 	${CMAKE_CURRENT_SOURCE_DIR}/table_printer.cpp
     )

##removing optimizations
SET_SOURCE_FILES_PROPERTIES(${BPRINTER_SOURCES} PROPERTIES COMPILE_FLAGS -O2)

###############################################################
add_definitions( -DBPRINTERLIBRARY )
add_definitions( -w )
#perform platform specific defintions
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DCPU86 )
    add_definitions( -DNO_TIMER )
else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DLINUX )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

###############################################################
add_library(bprinter STATIC ${BPRINTER_SOURCES})
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
target_link_libraries(bprinter m)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# install(TARGETS bprinter DESTINATION libs )
