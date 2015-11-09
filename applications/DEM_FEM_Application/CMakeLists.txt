set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDEM_FEM_Application ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
## include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../DEM_application )


## generate variables with the sources
set( KRATOS_DEM_FEM_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/dem_fem__application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dem_fem__python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosDEM_FEM_Application SHARED ${KRATOS_DEM_FEM_APPLICATION_SOURCES})
target_link_libraries(KratosDEM_FEM_Application KratosCore )
set_target_properties(KratosDEM_FEM_Application PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDEM_FEM_Application PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDEM_FEM_Application PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosDEM_FEM_Application DESTINATION libs )


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message("KratosDEM_FEM_Application subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DEM_FEM_Application.py" DESTINATION KratosMultiphysics )
