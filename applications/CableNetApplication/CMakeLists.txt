set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCableNetApplication ****")

################### PYBIND11
include(pybind11Tools)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## CableNet Core sources
file(GLOB_RECURSE KRATOS_CABLE_NET_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/cable_net_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/cable_net_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## CableNet testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CABLE_NET_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## CableNet python interface sources
file(GLOB_RECURSE KRATOS_CABLE_NET_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosCableNetCore SHARED ${KRATOS_CABLE_NET_APPLICATION_CORE} ${KRATOS_CABLE_NET_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosCableNetCore PRIVATE KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosCableNetCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_CABLE_NET_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosCableNetApplication MODULE THIN_LTO ${KRATOS_CABLE_NET_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosCableNetApplication PRIVATE KratosCableNetCore)
set_target_properties(KratosCableNetApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCableNetApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosCableNetApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosCableNetCore)
    cotire(KratosCableNetApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/CableNetApplication.py KratosMultiphysics/CableNetApplication/__init__.py )

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
install(TARGETS KratosCableNetCore DESTINATION libs )
install(TARGETS KratosCableNetApplication DESTINATION libs )
