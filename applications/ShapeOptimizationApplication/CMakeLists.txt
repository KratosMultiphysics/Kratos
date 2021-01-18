set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosShapeOptimizationApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## ShapeOptimization core sources
set( KRATOS_SHAPE_OPTIMIZATION_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/shape_optimization_application.cpp
)

## ShapeOptimization python interface sources
set( KRATOS_SHAPE_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/shape_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

add_library(KratosShapeOptimizationCore SHARED ${KRATOS_SHAPE_OPTIMIZATION_APPLICATION_CORE})
target_link_libraries(KratosShapeOptimizationCore PUBLIC KratosCore)
set_target_properties(KratosShapeOptimizationCore PROPERTIES COMPILE_DEFINITIONS "SHAPE_OPTIMIZATION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosShapeOptimizationApplication MODULE THIN_LTO ${KRATOS_SHAPE_OPTIMIZATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosShapeOptimizationApplication PRIVATE KratosShapeOptimizationCore)
set_target_properties(KratosShapeOptimizationApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosShapeOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosShapeOptimizationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosShapeOptimizationCore)
    cotire(KratosShapeOptimizationApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ShapeOptimizationApplication.py KratosMultiphysics/ShapeOptimizationApplication/__init__.py )

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
install(TARGETS KratosShapeOptimizationCore DESTINATION libs )
install(TARGETS KratosShapeOptimizationApplication DESTINATION libs )
