set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosTopologyOptimizationApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )

## generate variables with the sources
set( KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/topology_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_simp_element.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosTopologyOptimizationApplication SHARED ${KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_SOURCES})
target_link_libraries(KratosTopologyOptimizationApplication KratosCore KratosSolidMechanicsApplication)
set_target_properties(KratosTopologyOptimizationApplication PROPERTIES PREFIX "")
install(TARGETS KratosTopologyOptimizationApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosTopologyOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TopologyOptimizationApplication.py" DESTINATION KratosMultiphysics )
