set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMPMApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

# Define application core sources (pure C++ parts)
file(GLOB_RECURSE KRATOS_MPM_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/mpm_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mpm_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
)

# If there is mpi support, use..
if((${USE_MPI} MATCHES ON) )
  file( GLOB KRATOS_MPM_APPLICATION_CORE ${KRATOS_MPM_APPLICATION_CORE}
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mpi/*.cpp
  )
endif((${USE_MPI} MATCHES ON) )

# Define python interface sources
file(GLOB_RECURSE KRATOS_MPM_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

# Generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE  KRATOS_MPM_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMPMCore SHARED ${KRATOS_MPM_APPLICATION_CORE} ${KRATOS_MPM_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMPMCore PUBLIC KratosCore)
set_target_properties(KratosMPMCore PROPERTIES COMPILE_DEFINITIONS "MPM_APPLICATION=EXPORT,API")

# Define Kratos library which defines the basic python interface
pybind11_add_module(KratosMPMApplication MODULE THIN_LTO ${KRATOS_MPM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMPMApplication PRIVATE KratosMPMCore)
set_target_properties(KratosMPMApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMPMCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMPMApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMPMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMPMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MPMApplication.py KratosMultiphysics/MPMApplication/__init__.py )

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
install(TARGETS KratosMPMCore DESTINATION libs )
install(TARGETS KratosMPMApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMPMCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMPMApplication" PARENT_SCOPE)
