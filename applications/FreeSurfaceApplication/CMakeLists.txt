set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFreeSurfaceApplication ****")

include(pybind11Tools)
include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_FREE_SURFACE_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/free_surface_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/free_surface_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/free_surface_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_edgebased_levelset_solver_to_python.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosFreeSurfaceApplication MODULE ${KRATOS_FREE_SURFACE_APPLICATION_SOURCES})
target_link_libraries(KratosFreeSurfaceApplication KratosCore )
set_target_properties(KratosFreeSurfaceApplication PROPERTIES PREFIX "")
install(TARGETS KratosFreeSurfaceApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFreeSurfaceApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFreeSurfaceApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(USE_COTIRE MATCHES ON)
    cotire(KratosFreeSurfaceApplication)
endif(USE_COTIRE MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FreeSurfaceApplication.py" DESTINATION KratosMultiphysics )
