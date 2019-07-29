set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE
    ## MAIN FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application_variables.cpp
)

file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/*.cpp)
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE
        KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/bits/*.cpp
    )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## python interface sources
file(GLOB_RECURSE KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosStructuralMechanicsCore SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosStructuralMechanicsCore PUBLIC KratosCore)
set_target_properties(KratosStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosStructuralMechanicsApplication PUBLIC KratosStructuralMechanicsCore)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

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

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosStructuralMechanicsCore)
    cotire(KratosStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py" DESTINATION "KratosMultiphysics/StructuralMechanicsApplication" RENAME "__init__.py")
