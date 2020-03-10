set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosGeodataProcessingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## GeodataProcessing Core sources
file(GLOB_RECURSE KRATOS_GEODATA_PROCESSING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/geodata_processing_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geodata_processing_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## GeodataProcessing python interface sources
file(GLOB_RECURSE KRATOS_GEODATA_PROCESSING_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# GeodataProcessing testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_GEODATA_PROCESSING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosGeodataProcessingCore SHARED ${KRATOS_GEODATA_PROCESSING_APPLICATION_CORE} ${KRATOS_GEODATA_PROCESSING_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosGeodataProcessingCore PUBLIC KratosCore)
set_target_properties(KratosGeodataProcessingCore PROPERTIES COMPILE_DEFINITIONS "GEODATA_PROCESSING_APPLICATION=EXPORT,API")

# Define library Kratos which defines the basic python interface
pybind11_add_module(KratosGeodataProcessingApplication MODULE THIN_LTO ${KRATOS_GEODATA_PROCESSING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosGeodataProcessingApplication PUBLIC KratosGeodataProcessingCore)
set_target_properties(KratosGeodataProcessingApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosGeodataProcessingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosGeodataProcessingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosGeodataProcessingCore)
    cotire(KratosGeodataProcessingApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/GeodataProcessingApplication.py KratosMultiphysics/GeodataProcessingApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
)

# Install targets
install(TARGETS KratosGeodataProcessingCore DESTINATION libs )
install(TARGETS KratosGeodataProcessingApplication DESTINATION libs )