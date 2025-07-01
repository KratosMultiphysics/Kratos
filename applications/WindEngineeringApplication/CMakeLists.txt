set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosWindEngineeringApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## WindEngineering Core sources
file(GLOB_RECURSE KRATOS_WIND_ENGINEERING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/wind_engineering_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/wind_engineering_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## WindEngineering python interface sources
file(GLOB_RECURSE KRATOS_WIND_ENGINEERING_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosWindEngineeringCore SHARED ${KRATOS_WIND_ENGINEERING_APPLICATION_CORE})
target_link_libraries(KratosWindEngineeringCore PUBLIC KratosCore)
set_target_properties(KratosWindEngineeringCore PROPERTIES COMPILE_DEFINITIONS "WIND_ENGINEERING_APPLICATION=EXPORT,API")

## WindEngineering testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_WIND_ENGINEERING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
    kratos_add_gtests(TARGET KratosWindEngineeringCore SOURCES "${KRATOS_WIND_ENGINEERING_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosWindEngineeringApplication MODULE THIN_LTO ${KRATOS_WIND_ENGINEERING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosWindEngineeringApplication PRIVATE KratosWindEngineeringCore)
set_target_properties(KratosWindEngineeringApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosWindEngineeringApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosWindEngineeringApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/WindEngineeringApplication.py KratosMultiphysics/WindEngineeringApplication/__init__.py )

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
install(TARGETS KratosWindEngineeringCore DESTINATION libs )
install(TARGETS KratosWindEngineeringApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosWindEngineeringCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosWindEngineeringApplication" PARENT_SCOPE)