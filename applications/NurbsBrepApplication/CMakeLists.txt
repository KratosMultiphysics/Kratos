set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNurbBrepApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication ) 
#include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication ) 
#include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

# generate variables with the sources
set( KRATOS_NURBS_BREP_APPLICATION_SOURCES
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/nurbs_brep_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/nurbs_brep_application_variables.cpp
  
  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/nurbs_brep_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  
  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/FaceTrim.h
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/Edge.h
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/TrimmingCurve.h
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/BrepModel.h
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/Face.h  
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/BrepModelGeometryReader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/ControlPoint.h    
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/NurbsBrepModeler.cpp
    
  ##ELEMENTS
  # Adding the meshless base element
  #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_base_element.cpp
    
  ##CONDITIONS  ContinuityConditionPenalty
  #${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_base_condition.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosNurbsBrepApplication SHARED ${KRATOS_NURBS_BREP_APPLICATION_SOURCES})
target_link_libraries(KratosNurbsBrepApplication KratosStructuralMechanicsApplication KratosCore )
set_target_properties(KratosNurbsBrepApplication PROPERTIES PREFIX "")
install(TARGETS KratosNurbsBrepApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosNurbsBrepApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(KratosNURBSBRepApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/NURBSBRepApplication.py" DESTINATION KratosMultiphysics )
