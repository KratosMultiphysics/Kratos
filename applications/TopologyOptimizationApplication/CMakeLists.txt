set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosTopologyOptimizationApplication ****")

################### PYBIND11
include(pybind11Tools)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/TopologyOptimizationApplication )

## TopologyOptimization core sources
set( KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_CORE     
	${CMAKE_CURRENT_SOURCE_DIR}/topology_optimization_application.cpp 
    	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_simp_element.cpp
	)

## TopoloyOptimization python interface sources
set( KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE	
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/topology_optimization_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_simp_element.cpp
)

add_library(KratosTopologyOptimizationCore SHARED ${KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_CORE})
target_link_libraries(KratosTopologyOptimizationCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosTopologyOptimizationCore PROPERTIES COMPILE_DEFINITIONS "TOPOLOGY_OPTIMIZATION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosTopologyOptimizationApplication MODULE THIN_LTO ${KRATOS_TOPOLOGY_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosTopologyOptimizationApplication PRIVATE KratosTopologyOptimizationCore)
set_target_properties(KratosTopologyOptimizationApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosTopologyOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosTopologyOptimizationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/TopologyOptimizationApplication.py KratosMultiphysics/TopologyOptimizationApplication/__init__.py )

# Install python files neu dazu !!!
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosTopologyOptimizationCore DESTINATION libs )
install(TARGETS KratosTopologyOptimizationApplication DESTINATION libs )
