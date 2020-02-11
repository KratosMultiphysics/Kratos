set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMORApplication ****")

################### PYBIND11
include(pybind11Tools)

if( NOT DEFINED EIGEN_ROOT )
    message( FATAL_ERROR "EIGEN_ROOT not defined" )
endif()

add_definitions( -DEIGEN_DEFAULT_TO_ROW_MAJOR -DEIGEN_MPL2_ONLY )

include_directories(
    ${CMAKE_SOURCE_DIR}/kratos
        SYSTEM ${EIGEN_ROOT} # Including as System to suppress compile-warnings from it
)

# generate variables with the sources
set( KRATOS_MOR_APPLICATION_SOURCES
    ## MAIN FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/mor_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mor_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/acoustic_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/displacement_output_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/generalized_eigenvalue_utility.cpp

    ## PYTHON
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/mor_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosMORApplication MODULE THIN_LTO ${KRATOS_MOR_APPLICATION_SOURCES})
target_link_libraries(KratosMORApplication PRIVATE KratosCore )
set_target_properties(KratosMORApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosMORApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMORApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMORApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMORApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

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

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MORApplication.py" DESTINATION "KratosMultiphysics/MORApplication" RENAME "__init__.py")
