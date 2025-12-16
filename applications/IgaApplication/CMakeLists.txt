set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIgaApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories(
    ${KRATOS_SOURCE_DIR}/kratos
)

# Define application core sources (pure C++ parts)
file(GLOB_RECURSE KRATOS_IGA_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/iga_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
)

# Define python interface sources
file(GLOB_RECURSE KRATOS_IGA_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# Generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_IGA_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
  kratos_add_gtests(TARGET KratosIgaCore SOURCES "${KRATOS_IGA_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosIgaCore SHARED ${KRATOS_IGA_APPLICATION_CORE})
target_link_libraries(KratosIgaCore PUBLIC KratosCore)
set_target_properties(KratosIgaCore PROPERTIES COMPILE_DEFINITIONS "IGA_APPLICATION=EXPORT,API")

# Define library Kratos which defines the basic python interface
pybind11_add_module(KratosIgaApplication MODULE THIN_LTO ${KRATOS_IGA_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosIgaApplication PRIVATE KratosIgaCore)
set_target_properties(KratosIgaApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosIgaCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosIgaApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosIgaApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/IgaApplication.py KratosMultiphysics/IgaApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosIgaCore DESTINATION libs )
install(TARGETS KratosIgaApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosIgaCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosIgaApplication" PARENT_SCOPE)