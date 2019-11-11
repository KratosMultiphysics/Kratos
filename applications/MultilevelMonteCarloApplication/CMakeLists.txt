set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultilevelMonteCarloApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

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

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosMultilevelMonteCarloCore)
    cotire(KratosMultilevelMonteCarloApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries
install(TARGETS KratosMultilevelMonteCarloCore DESTINATION libs )
install(TARGETS KratosMultilevelMonteCarloApplication DESTINATION libs )

###############################################################################
# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultilevelMonteCarloApplication.py" DESTINATION "KratosMultiphysics/MultilevelMonteCarloApplication" RENAME "__init__.py")

###############################################################################
if(NOT DEFINED USING_PYCOMPSS)
    install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/exaqute" DESTINATION .)
endif(NOT DEFINED USING_PYCOMPSS)