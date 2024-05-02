set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDigitalTwinApplication ****")

find_package(FLANN REQUIRED)

################### PYBIND11
include(pybind11Tools)
include_directories( ${FLANN_INCLUDE_DIRS} )
include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## DigitalTwin Core sources
file(GLOB_RECURSE KRATOS_DIGITAL_TWIN_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/digital_twin_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/digital_twin_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_sensors/*.cpp

)

## DigitalTwin testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_DIGITAL_TWIN_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## DigitalTwin python interface sources
file(GLOB_RECURSE KRATOS_DIGITAL_TWIN_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

find_library(LZ4_LIBRARIES NAMES liblz4.so HINTS "lz4/lib")

add_library(KratosDigitalTwinCore SHARED ${KRATOS_DIGITAL_TWIN_APPLICATION_CORE} ${KRATOS_DIGITAL_TWIN_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosDigitalTwinCore PUBLIC KratosCore PRIVATE ${FLANN_LIBRARIES} ${LZ4_LIBRARIES})
set_target_properties(KratosDigitalTwinCore PROPERTIES COMPILE_DEFINITIONS "DIGITAL_TWIN_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDigitalTwinApplication MODULE THIN_LTO ${KRATOS_DIGITAL_TWIN_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDigitalTwinApplication PRIVATE KratosDigitalTwinCore)
set_target_properties(KratosDigitalTwinApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDigitalTwinApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosDigitalTwinApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/DigitalTwinApplication.py KratosMultiphysics/DigitalTwinApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/DigitalTwinApplication/python_registry_lists.py )

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
install(TARGETS KratosDigitalTwinCore DESTINATION libs )
install(TARGETS KratosDigitalTwinApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosDigitalTwinCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosDigitalTwinApplication" PARENT_SCOPE)