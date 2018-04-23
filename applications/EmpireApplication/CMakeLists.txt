set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosEmpireApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )


## generate variables with the sources
set( KRATOS_EMPIRE_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/empire_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/empire_python_application.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosEmpireApplication MODULE THIN_LTO ${KRATOS_EMPIRE_APPLICATION_SOURCES})
target_link_libraries(KratosEmpireApplication PRIVATE KratosCore )
set_target_properties(KratosEmpireApplication PROPERTIES PREFIX "")
install(TARGETS KratosEmpireApplication DESTINATION libs )


# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosEmpireApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosEmpireApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/EmpireApplication.py" DESTINATION KratosMultiphysics )
