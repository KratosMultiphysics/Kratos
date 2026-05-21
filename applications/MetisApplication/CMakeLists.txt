set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMetisApplication ****")

include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

find_package(METIS REQUIRED)
include_directories(${METIS_INCLUDE_DIR})
link_directories(${METIS_LIBRARIES})
install(FILES ${METIS_LIBRARIES} DESTINATION libs)

## Metis Core sources
file( GLOB KRATOS_METIS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/metis_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## Metis python interface sources
file(GLOB_RECURSE KRATOS_METIS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## Metis Testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB KRATOS_METIS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(TARGET KratosMetisCore SOURCES "${KRATOS_METIS_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMetisCore SHARED ${KRATOS_METIS_APPLICATION_CORE})
target_link_libraries(KratosMetisCore PUBLIC KratosCore ${METIS_LIBRARIES})
set_target_properties(KratosMetisCore PROPERTIES COMPILE_DEFINITIONS "METIS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMetisApplication MODULE THIN_LTO ${KRATOS_METIS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMetisApplication PRIVATE KratosMetisCore)
set_target_properties(KratosMetisApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMetisCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMetisApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMetisApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMetisApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MetisApplication.py KratosMultiphysics/MetisApplication/__init__.py )

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
install(TARGETS KratosMetisCore DESTINATION libs )
install(TARGETS KratosMetisApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMetisCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMetisApplication" PARENT_SCOPE)
