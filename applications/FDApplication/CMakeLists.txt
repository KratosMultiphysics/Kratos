set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFDApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_FD_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/fd_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fd_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/grid_printer.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/bfecc.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFDApplication SHARED ${KRATOS_FD_APPLICATION_SOURCES})
target_link_libraries(KratosFDApplication KratosCore )
set_target_properties(KratosFDApplication PROPERTIES PREFIX "")
install(TARGETS KratosFDApplication DESTINATION libs )

if(${ACTIVATE_DEBUG_MACRO} MATCHES ON)        #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosFDApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFDApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FDApplication.py" DESTINATION KratosMultiphysics )
