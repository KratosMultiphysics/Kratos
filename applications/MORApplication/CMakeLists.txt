set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMORApplication ****")

################### PYBIND11
include(pybind11Tools)

## add Eigen sources from the LinearSolversApplication
add_definitions( -DEIGEN_DEFAULT_TO_ROW_MAJOR -DEIGEN_MPL2_ONLY )
include_directories(
    ${CMAKE_SOURCE_DIR}/kratos
        SYSTEM ../LinearSolversApplication/external_libraries/eigen3 # Including as System to suppress compile-warnings from it
)

## MOR core sources
file(GLOB_RECURSE KRATOS_MOR_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/mor_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mor_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## MOR python interface sources
file(GLOB_RECURSE KRATOS_MOR_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MOR_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMORCore SHARED ${KRATOS_MOR_APPLICATION_CORE} ${KRATOS_MOR_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMORCore PUBLIC KratosCore)
set_target_properties(KratosMORCore PROPERTIES COMPILE_DEFINITIONS "MOR_APPLICATION=EXPORT,API")

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMORApplication MODULE THIN_LTO ${KRATOS_MOR_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMORApplication PUBLIC KratosMORCore )
set_target_properties(KratosMORApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosMORApplication)
endif(USE_COTIRE MATCHES ON)


# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMORApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMORApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MORApplication.py KratosMultiphysics/MORApplication/__init__.py )

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

# Add to the KratosMultiphisics Python module
install(TARGETS KratosMORCore DESTINATION libs )
install(TARGETS KratosMORApplication DESTINATION libs )
