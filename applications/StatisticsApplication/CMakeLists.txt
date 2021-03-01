set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStatisticsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Statistics Core sources
file(GLOB_RECURSE KRATOS_STATISTICS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/statistics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/statistics_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_methods/*.cpp
)

## Statistics testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
file(GLOB_RECURSE KRATOS_STATISTICS_APPLICATION_TESTING_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp
)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Statistics python interface sources
file(GLOB_RECURSE KRATOS_STATISTICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosStatisticsCore SHARED ${KRATOS_STATISTICS_APPLICATION_CORE} ${KRATOS_STATISTICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosStatisticsCore PUBLIC KratosCore)
set_target_properties(KratosStatisticsCore PROPERTIES COMPILE_DEFINITIONS "STATISTICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosStatisticsApplication MODULE THIN_LTO ${KRATOS_STATISTICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStatisticsApplication PRIVATE KratosStatisticsCore)
set_target_properties(KratosStatisticsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosStatisticsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosStatisticsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosStatisticsCore)
    cotire(KratosStatisticsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/StatisticsApplication.py KratosMultiphysics/StatisticsApplication/__init__.py )

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
install(TARGETS KratosStatisticsCore DESTINATION libs )
install(TARGETS KratosStatisticsApplication DESTINATION libs )
