set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosShallowWaterApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
file( GLOB_RECURSE KRATOS_SHALLOW_WATER_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/shallow_water_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/shallow_water_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## define python interface sources
file( GLOB_RECURSE KRATOS_SHALLOW_WATER_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## define shallow water test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file( GLOB_RECURSE KRATOS_SHALLOW_WATER_TESTING
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
    )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## ShallowWaterApplication core library (C++ parts)
add_library( KratosShallowWaterCore SHARED ${KRATOS_SHALLOW_WATER_APPLICATION_SOURCES} ${KRATOS_SHALLOW_WATER_TESTING})
target_link_libraries( KratosShallowWaterCore PUBLIC KratosCore)
set_target_properties( KratosShallowWaterCore PROPERTIES COMPILE_DEFINITIONS "SHALLOW_WATER_APPLICATION=EXPORT,API")

## define library Kratos which defines the basic python interface
pybind11_add_module( KratosShallowWaterApplication MODULE THIN_LTO ${KRATOS_SHALLOW_WATER_APPLICATION_PYTHON_INTERFACE})
target_link_libraries( KratosShallowWaterApplication PRIVATE KratosShallowWaterCore)
set_target_properties( KratosShallowWaterApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosShallowWaterApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosShallowWaterApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosShallowWaterCore)
    cotire(KratosShallowWaterApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ShallowWaterApplication.py KratosMultiphysics/ShallowWaterApplication/__init__.py )

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
install( TARGETS KratosShallowWaterCore DESTINATION libs )
install( TARGETS KratosShallowWaterApplication DESTINATION libs )
