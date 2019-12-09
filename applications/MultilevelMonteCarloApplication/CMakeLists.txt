set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultilevelMonteCarloApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/multilevel_monte_carlo_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multilevel_monte_carlo_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multilevel_monte_carlo_python_application.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMultilevelMonteCarloApplication MODULE THIN_LTO ${KRATOS_MULTILEVEL_MONTE_CARLO_APPLICATION_SOURCES})
target_link_libraries(KratosMultilevelMonteCarloApplication PRIVATE KratosCore )
set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMultilevelMonteCarloApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosMultilevelMonteCarloApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MultilevelMonteCarloApplication.py KratosMultiphysics/MultilevelMonteCarloApplication/__init__.py )

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
install(TARGETS KratosMultilevelMonteCarloApplication DESTINATION libs )

# Install pycompss mockup lib
if(NOT DEFINED USING_PYCOMPSS)
    kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/exaqute KratosMultiphysics/${CURRENT_DIR_NAME} )
endif(NOT DEFINED USING_PYCOMPSS)