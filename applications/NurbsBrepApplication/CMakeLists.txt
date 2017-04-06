set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNurbBrepApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

file(GLOB KRATOS_NURBS_BREP_APPLICATION_PYTHON_SOURCES 
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp")
file(GLOB KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES 
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.hpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/knot_span/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep/*.cpp")
file(GLOB KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

#groups
#source_group("Main" FILES ${KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES})
#source_group("Python" FILES ${KRATOS_NURBS_BREP_APPLICATION_PYTHON_SOURCES})
#source_group("Brep" FILES ${KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES})


# define library Kratos which defines the basic python interface
add_library(KratosNurbsBrepApplication SHARED ${KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES}
                                              ${KRATOS_NURBS_BREP_APPLICATION_PYTHON_SOURCES}
                                              ${KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES})
target_link_libraries(KratosNurbsBrepApplication KratosCore )
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
