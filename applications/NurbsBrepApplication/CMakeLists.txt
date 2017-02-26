set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNurbBrepApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
file(GLOB KRATOS_NURBS_BREP_APPLICATION_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp")

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
