set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring KratosContactStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ContactStructuralMechanicsApplication )

## ContactStructuralMechanics Core sources
file(GLOB_RECURSE KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_CORE
    ## MAIN FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application_variables.cpp
    ## UTILITIES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
    ## FRICTIONAL LAWS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_frictional_laws/*.cpp
    ## PROCESSES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ## CONDITIONS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ## CONSTRAINTS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_master_slave_constraints/*.cpp
    ## EXTERNAL LIBRARIES
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ContactStructuralMechanics python interface sources
file(GLOB_RECURSE KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosContactStructuralMechanicsCore SHARED ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosContactStructuralMechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosContactStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "CONTACT_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosContactStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosContactStructuralMechanicsApplication PRIVATE KratosContactStructuralMechanicsCore)
set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    ## The interface_preprocess is not compatible with cotire (solve problem with explicit instantations)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_preprocess.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    cotire(KratosContactStructuralMechanicsCore)
    cotire(KratosContactStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ContactStructuralMechanicsApplication.py KratosMultiphysics/ContactStructuralMechanicsApplication/__init__.py )

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
install(TARGETS KratosContactStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosContactStructuralMechanicsApplication DESTINATION libs )