set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosLaserDrillingApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_LASERDRILLING_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/laserdrilling_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/laserdrilling_application_variables.cpp
)

## define python interface sources
set( KRATOS_LASERDRILLING_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_laserdrilling_python_application.cpp
)

## define test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
     file(GLOB_RECURSE KRATOS_LASERDRILLING_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## LaserDrillingApplication core library (C++ parts)
add_library( KratosLaserDrillingCore SHARED ${KRATOS_LASERDRILLING_APPLICATION_CORE_SOURCES} ${KRATOS_LASERDRILLING_TESTING_SOURCES} )
target_link_libraries( KratosLaserDrillingCore PUBLIC KratosCore)
set_target_properties( KratosLaserDrillingCore PROPERTIES COMPILE_DEFINITIONS "LASERDRILLING_APPLICATION=EXPORT,API")

## LaserDrillingApplication python module
pybind11_add_module( KratosLaserDrillingApplication MODULE THIN_LTO ${KRATOS_LASERDRILLING_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosLaserDrillingApplication PRIVATE KratosLaserDrillingCore)
set_target_properties( KratosLaserDrillingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosLaserDrillingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosLaserDrillingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosLaserDrillingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosLaserDrillingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/LaserDrillingApplication.py KratosMultiphysics/LaserDrillingApplication/__init__.py )

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

# Install Targets
install(TARGETS KratosLaserDrillingCore DESTINATION libs )
install(TARGETS KratosLaserDrillingApplication DESTINATION libs )
