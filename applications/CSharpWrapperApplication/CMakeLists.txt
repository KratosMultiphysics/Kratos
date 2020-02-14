set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring CSharpWrapperApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories(${CMAKE_SOURCE_DIR}/kratos)
include_directories(${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication)

set(KRATOS_WRAPPER_SOURCES_CORE
        ## MAIN FILES
        ${CMAKE_CURRENT_SOURCE_DIR}/csharp_wrapper_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/csharp_wrapper_application_variables.cpp

        ## WRAPPER API
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/kratos_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/mesh_converter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/vector3.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/kratos_internals.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/id_translator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/kratos_wrapper.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_sources/model_part_wrapper.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/model_part_wrapper_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/node_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/element_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/condition_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/variable_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/utils_external_bindings.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/external_bindings/id_translator_external_bindings.cpp)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_WRAPPER_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## generate variables with the sources
set( KRATOS_WRAPPER_SOURCES_PYTHON_INTERFACE
        ## CUSTOM PYTHON
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/csharp_wrapper_python_application.cpp
        )

add_library(KratosCSharpWrapperCore SHARED ${KRATOS_WRAPPER_SOURCES_CORE} ${KRATOS_WRAPPER_TESTING_SOURCES})
target_link_libraries(KratosCSharpWrapperCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosCSharpWrapperCore PROPERTIES COMPILE_DEFINITIONS "CSHARP_WRAPPER_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosCSharpWrapperApplication MODULE THIN_LTO ${KRATOS_WRAPPER_SOURCES_PYTHON_INTERFACE})
target_link_libraries(KratosCSharpWrapperApplication PRIVATE KratosCSharpWrapperCore)
set_target_properties(KratosCSharpWrapperApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCSharpWrapperApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosCSharpWrapperApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosCSharpWrapperCore)
    cotire(KratosCSharpWrapperApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/CSharpWrapperApplication.py KratosMultiphysics/CSharpWrapperApplication/__init__.py )

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
install(TARGETS KratosCSharpWrapperCore DESTINATION libs)
install(TARGETS KratosCSharpWrapperApplication DESTINATION libs)