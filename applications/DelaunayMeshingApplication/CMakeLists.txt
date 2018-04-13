set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDelaunayMeshingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_DELAUNAY_MESHING_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/petsc_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/petsc_application_variables.cpp
)

## Generate variables with the sources
set( KRATOS_DELAUNAY_MESHING_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/petsc_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

add_library(KratosDelaunayMeshingCore SHARED ${KRATOS_DELAUNAY_MESHING_APPLICATION_CORE})
target_link_libraries(KratosDelaunayMeshingCore PUBLIC KratosCore)
set_target_properties(KratosDelaunayMeshingCore PROPERTIES COMPILE_DEFINITIONS "DELAUNAY_MESHING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDelaunayMeshingApplication MODULE ${KRATOS_DELAUNAY_MESHING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDelaunayMeshingApplication PUBLIC KratosDelaunayMeshingCore)
set_target_properties(KratosDelaunayMeshingApplication PROPERTIES PREFIX "")


# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDelaunayMeshingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
set_target_properties(KratosDelaunayMeshingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install testing files
if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    cotire(KratosDelaunayMeshingCore)
    cotire(KratosDelaunayMeshingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosDelaunayMeshingCore DESTINATION libs )
install(TARGETS KratosDelaunayMeshingApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DelaunayMeshingApplication.py" DESTINATION KratosMultiphysics )
