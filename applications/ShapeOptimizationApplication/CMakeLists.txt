set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosShapeOptimizationApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ExternalSolversApplication )

if(${INCLUDE_FEAST} MATCHES ON)
#  include(CMakeAddFortranSubdirectory)
#  cmake_add_fortran_subdirectory(${CMAKE_SOURCE_DIR}/applications/ExternalSolversApplication/custom_external_libraries/FEAST NO_EXTERNAL_INSTALL )
#  if (NOT DEFINED FEAST_ROOT)
#    set(FEAST_ROOT ${CMAKE_SOURCE_DIR}/applications/ExternalSolversApplication/custom_external_libraries/FEAST/3.0 )
#  endif(NOT DEFINED FEAST_ROOT)
  include_directories( ${CMAKE_SOURCE_DIR}/applications/ExternalSolversApplication/custom_external_libraries/FEAST/3.0/include )
  add_definitions(-DINCLUDE_FEAST)
endif(${INCLUDE_FEAST} MATCHES ON)

## generate variables with the sources
set( KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/shape_optimization_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/shape_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/shape_optimization_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_analytic_sensitivity_element.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosShapeOptimizationApplication SHARED ${KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES})
target_link_libraries(KratosShapeOptimizationApplication KratosCore KratosSolidMechanicsApplication KratosExternalSolversApplication)
set_target_properties(KratosShapeOptimizationApplication PROPERTIES PREFIX "")
install(TARGETS KratosShapeOptimizationApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosShapeOptimizationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ShapeOptimizationApplication.py" DESTINATION KratosMultiphysics )
