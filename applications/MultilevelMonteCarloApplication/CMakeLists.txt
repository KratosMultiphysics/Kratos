set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultilevelMonteCarloApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
set( KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_SOURCES
    # application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/multilevel_monte_carlo_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multilevel_monte_carlo_application_variables.cpp
    # statistics
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_statistics/power_sums_statistics.cpp
)

## define python interface sources
set (KRATOS_MULTILEVEL_MONTE_CARLO_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multilevel_monte_carlo_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_statistics_to_python.cpp
)

## define multilevel monte carlo test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## MultilevelMonteCarloApplication core library (C++ parts)
add_library(KratosMultilevelMonteCarloCore SHARED ${KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_SOURCES} ${KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMultilevelMonteCarloCore PUBLIC KratosCore)
set_target_properties(KratosMultilevelMonteCarloCore PROPERTIES COMPILE_DEFINITIONS "MULTILEVEL_MONTE_CARLO_APPLICATION=EXPORT,API")

## MultilevelMonteCarloApplication python module
pybind11_add_module(KratosMultilevelMonteCarloApplication MODULE THIN_LTO ${KRATOS_MULTILEVEL_MONTE_CARLO_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosMultilevelMonteCarloApplication PRIVATE KratosMultilevelMonteCarloCore)
set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMultilevelMonteCarloCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

## changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

## changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

## add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MultilevelMonteCarloApplication.py KratosMultiphysics/MultilevelMonteCarloApplication/__init__.py )

## install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

## Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
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

## install targets
install(TARGETS KratosMultilevelMonteCarloCore DESTINATION libs )
install(TARGETS KratosMultilevelMonteCarloApplication DESTINATION libs )

############################################################################
#####
#####        PyCOMPSs (Distributed Computation library)
#####
############################################################################
# Install PyCOMPSs mockup library
if(NOT DEFINED USING_PYCOMPSS OR NOT USING_PYCOMPSS)
    message ( STATUS "Automatic installation of PyCOMPSs mockup library" )
    add_definitions(-DKRATOS_MLMC_USE_MOCKUP_PYCOMPSS=1)
	kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/PyCOMPSs KratosMultiphysics/${CURRENT_DIR_NAME}/PyCOMPSs )
else (USING_PYCOMPSS)
    message ( STATUS "Manual installation of PyCOMPSs library" )
    add_definitions(-DKRATOS_MLMC_USE_MOCKUP_PYCOMPSS=0)
endif(NOT DEFINED USING_PYCOMPSS OR NOT USING_PYCOMPSS)

############################################################################
#####
#####        XMC (Uncertainty Quantification library)
#####
############################################################################
# Install XMC library
if (NOT DEFINED USING_INTERNAL_XMC OR USING_INTERNAL_XMC)
    message ( STATUS "Automatic installation of XMC library" )
    add_definitions(-DKRATOS_MLMC_USE_INTERNAL_XMC=1)
	kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/XMC KratosMultiphysics/${CURRENT_DIR_NAME}/XMC )
else (NOT USING_INTERNAL_XMC)
    message ( STATUS "Manual installation of XMC library" )
    add_definitions(-DKRATOS_MLMC_USE_INTERNAL_XMC=0)
endif(NOT DEFINED USING_INTERNAL_XMC OR USING_INTERNAL_XMC)

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMultilevelMonteCarloCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMultilevelMonteCarloApplication" PARENT_SCOPE)