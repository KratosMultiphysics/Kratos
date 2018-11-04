set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosShapeOptimizationApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/shape_optimization_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/shape_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/shape_optimization_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosShapeOptimizationApplication MODULE THIN_LTO ${KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES})
target_link_libraries(KratosShapeOptimizationApplication PRIVATE KratosCore)
set_target_properties(KratosShapeOptimizationApplication PROPERTIES PREFIX "")
install(TARGETS KratosShapeOptimizationApplication DESTINATION libs )

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosShapeOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

###############################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosShapeOptimizationApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(${INSTALL_TESTING_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ShapeOptimizationApplication.py" DESTINATION KratosMultiphysics )