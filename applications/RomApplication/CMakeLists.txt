set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRomApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## RomApplication Core sources
file(GLOB_RECURSE KRATOS_ROM_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/rom_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rom_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## RomApplication testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_ROM_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## RomApplication python interface sources
file(GLOB_RECURSE KRATOS_ROM_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosRomApplicationCore SHARED ${KRATOS_ROM_APPLICATION_CORE} ${KRATOS_ROM_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosRomApplicationCore PUBLIC KratosCore)
set_target_properties(KratosRomApplicationCore PROPERTIES COMPILE_DEFINITIONS "ROM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosRomApplication MODULE THIN_LTO ${KRATOS_ROM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosRomApplication PRIVATE KratosRomApplicationCore)
set_target_properties(KratosRomApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRomApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRomApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/RomApplication.py KratosMultiphysics/RomApplication/__init__.py )

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

if(USE_COTIRE MATCHES ON)
    cotire(KratosRomApplicationCore)
    cotire(KratosRomApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosRomApplicationCore DESTINATION libs )
install(TARGETS KratosRomApplication DESTINATION libs )
