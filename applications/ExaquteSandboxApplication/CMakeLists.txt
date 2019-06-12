set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExaquteSandboxApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/exaqute_sandbox_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/exaqute_sandbox_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp


# define library Kratos which defines the basic python interface
pybind11_add_module(KratosExaquteSandboxApplication MODULE THIN_LTO ${KRATOS_EXAQUTE_SANDBOX_APPLICATION_SOURCES})
target_link_libraries(KratosExaquteSandboxApplication PRIVATE KratosCore )
set_target_properties(KratosExaquteSandboxApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosExaquteSandboxApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosExaquteSandboxApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosExaquteSandboxApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ExaquteSandboxApplication.py" DESTINATION "KratosMultiphysics/ExaquteSandboxApplication" RENAME "__init__.py")
