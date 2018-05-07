set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshMovingApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_MESH_MOVING_APPLICATION_SOURCES
	      ${CMAKE_CURRENT_SOURCE_DIR}/mesh_moving_application.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/laplacian_meshmoving_element.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/structural_meshmoving_element.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	      ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mesh_moving_python_application.cpp
		  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/move_mesh_utilities.cpp
		  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/explicit_mesh_moving_utilities.cpp
)

## define mesh moving test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_MESH_MOVING_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMeshMovingApplication MODULE THIN_LTO ${KRATOS_MESH_MOVING_APPLICATION_SOURCES} ${KRATOS_MESH_MOVING_TESTING_SOURCES})
target_link_libraries(KratosMeshMovingApplication PRIVATE KratosCore )
set_target_properties(KratosMeshMovingApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMeshMovingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMeshMovingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosMeshMovingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMeshMovingApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MeshMovingApplication.py" DESTINATION KratosMultiphysics )
