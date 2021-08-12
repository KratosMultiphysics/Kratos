set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring WindEngineeringApplication ****")

include_directories("${KRATOS_SOURCE_DIR}/kratos")

# Placeholder for setting up pybind

# Application library
set(KRATOS_WIND_ENGINEERING_APPLICATION_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/wind_engineering_application.cpp")

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_WIND_ENGINEERING_APPLICATION_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tets/*.cpp)
endif()

add_library(KratosWindEngineeringApplicationCore SHARED
            ${KRATOS_WIND_ENGINEERING_APPLICATION_SOURCES}
            ${KRATOS_WIND_ENGINEERING_APPLICATION_TEST_SOURCES})
target_link_libraries(KratosWindEngineeringApplicationCore PRIVATE KratosCore)

# Python bindings
include(pybind11Tools)

set(KRATOS_WIND_ENGINEERING_PYTHON_INTERFACE_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/custom_python/wind_engineering_python_application.cpp")

pybind11_add_module(KratosWindEngineeringApplication MODULE THIN_LTO ${KRATOS_WIND_ENGINEERING_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosWindEngineeringApplication PRIVATE KratosWindEngineeringApplicationCore)
set_target_properties(KratosWindEngineeringApplication PROPERTIES PREFIX "")

# Set batch size in unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosWindEngineeringApplicationCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosWindEngineeringApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Set the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosWindEngineeringApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# Set the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosWindEngineeringApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/WindEngineeringApplication.py KratosMultiphysics/WindEngineeringApplication/__init__.py )

# Install python scripts
get_filename_component(WIND_ENGINEERING_APPLICATION_SOURCE_DIR_NAME
                       "${CMAKE_CURRENT_SOURCE_DIR}"
                       NAME)

kratos_python_install_directory("${INSTALL_PYTHON_USING_LINKS}"
                                "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts"
                                "KratosMultiphysics/${WIND_ENGINEERING_APPLICATION_SOURCE_DIR_NAME}")

                                
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

# Install binaries
install(TARGETS KratosWindEngineeringApplicationCore DESTINATION libs)
install(TARGETS KratosWindEngineeringApplication DESTINATION libs)