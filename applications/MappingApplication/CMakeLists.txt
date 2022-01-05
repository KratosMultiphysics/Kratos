set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMappingApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Mapping Core sources
file( GLOB KRATOS_MAPPING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mapping_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_mappers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_searching/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/*.cpp
)

## define test sources (cpp tests)
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB KRATOS_MAPPING_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Mapping python interface sources
file(GLOB_RECURSE KRATOS_MAPPING_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosMappingCore SHARED ${KRATOS_MAPPING_APPLICATION_CORE} ${KRATOS_MAPPING_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMappingCore PUBLIC KratosCore )
set_target_properties(KratosMappingCore PROPERTIES COMPILE_DEFINITIONS "MAPPING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMappingApplication MODULE THIN_LTO ${KRATOS_MAPPING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMappingApplication PRIVATE KratosMappingCore)
set_target_properties(KratosMappingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMappingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMappingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMappingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Unity build
if(CMAKE_UNITY_BUILD MATCHES ON)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mapping_matrix_utilities.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_vector_container.cpp PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
endif(CMAKE_UNITY_BUILD MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MappingApplication.py KratosMultiphysics/MappingApplication/__init__.py )

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
install(TARGETS KratosMappingCore DESTINATION libs )
install(TARGETS KratosMappingApplication DESTINATION libs )

## If MPI is used, then enable the MPIExtension
if(${USE_MPI})
  message("Enabling mpi extension for MappingApplication")
  add_subdirectory(mpi_extension)
endif(${USE_MPI})

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMappingCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMappingApplication" PARENT_SCOPE)