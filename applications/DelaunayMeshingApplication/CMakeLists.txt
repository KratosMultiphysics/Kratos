set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDelaunayMeshingApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

################### PYBIND11
include(pybind11Tools)

## Set a default value for triangle and tetgen in case the user does not set it
## to avoid problems with the define linkage block. By default we will compile the lib
if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_DIR "${CMAKE_SOURCE_DIR}/external_libraries/triangle")
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/tetgen")
endif(NOT DEFINED ${INCLUDE_TETGEN})

set(MESHER_LIBRARIES "")
if(${INCLUDE_TRIANGLE} MATCHES ON)
  include_directories(${TRIANGLE_DIR})
  list(APPEND MESHER_LIBRARIES triangle)
endif(${INCLUDE_TRIANGLE} MATCHES ON)


if(${INCLUDE_TETGEN} MATCHES ON)
  set(EXTERNAL_LIBRARIES_DIR "${CMAKE_CURRENT_SOURCE_DIR}/external_libraries")
  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/custom_cmake/Tetgen.cmake")
  include_directories(${TETGEN_DIR})
  add_definitions(-DINCLUDE_TETGEN)
  add_definitions(-DTETLIBRARY)
  add_definitions(-DTETGEN150)
  add_definitions(-w)
  #perform platform specific defintions
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DCPU86)
    add_definitions(-DNO_TIMER)
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions(-DLINUX)
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set(TETGEN_SOURCES ${TETGEN_DIR}/tetgen.cxx)
  set(UNOPTIMIZED_PREDICATES ${TETGEN_DIR}/predicates.cxx)
  SET_SOURCE_FILES_PROPERTIES(${UNOPTIMIZED_PREDICATES} PROPERTIES COMPILE_FLAGS -O0)
  SET_SOURCE_FILES_PROPERTIES(${TETGEN_SOURCES} PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS})
  add_library(tetgen STATIC ${TETGEN_SOURCES} ${UNOPTIMIZED_PREDICATES})
  list(APPEND MESHER_LIBRARIES tetgen)
endif(${INCLUDE_TETGEN} MATCHES ON)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  list(APPEND MESHER_LIBRARIES m)
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")

## Generate variables with the sources
set( KRATOS_DELAUNAY_MESHING_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/delaunay_meshing_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/delaunay_meshing_application_variables.cpp

	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesher_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_data_transfer_utilities.cpp

	##conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/composite_condition.cpp

	##meshers
	${CMAKE_CURRENT_SOURCE_DIR}/custom_meshers/mesher.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_meshers/triangular_mesh_2D_mesher.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_meshers/tetrahedral_mesh_3D_mesher.cpp

)

## Generate variables with the sources
set( KRATOS_DELAUNAY_MESHING_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/delaunay_meshing_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_meshers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp
)

add_library(KratosDelaunayMeshingCore SHARED ${KRATOS_DELAUNAY_MESHING_APPLICATION_CORE})
target_link_libraries(KratosDelaunayMeshingCore PUBLIC KratosCore ${MESHER_LIBRARIES})
set_target_properties(KratosDelaunayMeshingCore PROPERTIES COMPILE_DEFINITIONS "DELAUNAY_MESHING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDelaunayMeshingApplication MODULE ${KRATOS_DELAUNAY_MESHING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDelaunayMeshingApplication PRIVATE KratosDelaunayMeshingCore)
set_target_properties(KratosDelaunayMeshingApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDelaunayMeshingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDelaunayMeshingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosDelaunayMeshingCore)
    cotire(KratosDelaunayMeshingApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/DelaunayMeshingApplication.py KratosMultiphysics/DelaunayMeshingApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosDelaunayMeshingCore DESTINATION libs )
install(TARGETS KratosDelaunayMeshingApplication DESTINATION libs )