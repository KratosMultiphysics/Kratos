set(CMAKE_INCLUDE_CURRENT_DIR ON)

if ( ${EMPIRE_NURBS_VERTEX_MORPHING} MATCHES ON )
  message("**** configuring KratosShapeOptimizationApplication with EMPIRE NURBS VERTEX MORPHING ****")
else ( ${EMPIRE_NURBS_VERTEX_MORPHING} MATCHES ON )
  message("**** configuring KratosShapeOptimizationApplication ****")
endif ( ${EMPIRE_NURBS_VERTEX_MORPHING} MATCHES ON )

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set( KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/shape_optimization_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/shape_optimization_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/shape_optimization_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_analytic_sensitivity_element.cpp
)

if ( ${EMPIRE_NURBS_VERTEX_MORPHING} MATCHES ON )
  add_definitions(-DEMPIRE_NURBS_VERTEX_MORPHING)
  #empire header files
  include_directories( $ENV{EMPIRE_MAPPER_LIB_INC_ON_MACHINE} )
  find_package(OCE)
  include_directories( ${OCE_INCLUDE_DIRS} )
endif(${EMPIRE_NURBS_VERTEX_MORPHING} MATCHES ON)

find_package(HDF5 REQUIRED COMPONENTS CXX)
include_directories( ${HDF5_INCLUDE_DIRS} )

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosShapeOptimizationApplication SHARED ${KRATOS_SHAPEOPTIMIZATION_APPLICATION_SOURCES})
target_link_libraries(KratosShapeOptimizationApplication KratosCore KratosStructuralMechanicsApplication ${HDF5_CXX_LIBRARIES} ${OCE_LIBRARIES} $ENV{EMPIRE_MAPPER_LIBSO_ON_MACHINE} )
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