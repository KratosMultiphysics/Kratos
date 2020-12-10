set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosTopologyOptimizationApplication ****")

include(pybind11Tools)
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
#include_directories( ${CMAKE_SOURCE_DIR}/applications/ )

#dependecies on other applications
#include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructualMechanicsApplication )

## generate variables with the sources
set( KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_CORE     
	${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization_application.cpp 
	)
	
set( KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/topology_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp

	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_simp_element.cpp
)
######################## neu

#############################
###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosTopologyOptimizationApplication MODULE THIN_LTO ${KRATOS_TOPOLOGYOPTIMIZATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosTopologyOptimizationApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosTopologyOptimizationApplication PROPERTIES PREFIX "")
#install(TARGETS KratosTopologyOptimizationApplication DESTINATION libs )


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
#install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/TopologyOptimizationApplication.py" DESTINATION KratosMultiphysics )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/TopologyOptimizationApplication.py KratosMultiphysics/TopologyOptimizationApplication/__init__.py )

# Install python files neu dazu !!!
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosTopologyOptimizationApplication DESTINATION libs )
