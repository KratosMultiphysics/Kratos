set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNurbBrepApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories(
  ${CMAKE_SOURCE_DIR}/kratos
  ${CMAKE_SOURCE_DIR}/external_libraries/ANurbs/include
)

file(GLOB KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_search/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/knot_span/*.cpp"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/brep/*.cpp")
file(GLOB KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.h"
                                                "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

## Define python interface sources
set( KRATOS_NURBS_BREP_APPLICATION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/nurbs_brep_python_application.cpp
)
												
#groups
#source_group("Main" FILES ${KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES})
#source_group("Python" FILES ${KRATOS_NURBS_BREP_APPLICATION_PYTHON_SOURCES})
#source_group("Brep" FILES ${KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES})


# define library Kratos which defines the basic python interface
#add_library(KratosNurbsBrepCore SHARED ${KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES}
#                                       ${KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES})
#target_link_libraries(KratosNurbsBrepCore PRIVATE KratosCore )
#set_target_properties(KratosNurbsBrepCore PROPERTIES COMPILE_DEFINITIONS "NURBS_BREP_APPLICATION=EXPORT,API")


pybind11_add_module(KratosNurbsBrepApplication MODULE THIN_LTO ${KRATOS_NURBS_BREP_APPLICATION_PYTHON_INTERFACE_SOURCES}
                                               ${KRATOS_NURBS_BREP_APPLICATION_MAIN_SOURCES}
                                       ${KRATOS_NURBS_BREP_APPLICATION_BREP_SOURCES})
target_link_libraries(KratosNurbsBrepApplication PRIVATE KratosCore )
set_target_properties(KratosNurbsBrepApplication PROPERTIES COMPILE_DEFINITIONS "IGA_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosNurbsBrepApplication PROPERTIES PREFIX "")



# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  set_target_properties(KratosNurbsBrepApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(KratosNurbsBrepApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    #cotire(KratosNurbsBrepCore)
    cotire(KratosNurbsBrepApplication)
endif(USE_COTIRE MATCHES ON)

#install(TARGETS KratosNurbsBrepCore DESTINATION libs )
install(TARGETS KratosNurbsBrepApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/NurbsBrepApplication.py" DESTINATION KratosMultiphysics )
