set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosViscosityModulatorApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_VISCOSITY_MODULATOR_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/viscosity_modulator_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/viscosity_modulator_application_variables.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/viscosity_modulator_element.cpp

    # Utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/viscosity_modulator_coupling_utilities.cpp

    # Process
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/boussinesq_modulator_field_process.cpp
)

## define python interface sources
set( KRATOS_VISCOSITY_MODULATOR_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_viscosity_modulator_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

## define test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file( GLOB_RECURSE KRATOS_VISCOSITY_MODULATOR_TESTING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
    )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## ViscosityModulatorApplication core library (C++ parts)
add_library( KratosViscosityModulatorCore SHARED ${KRATOS_VISCOSITY_MODULATOR_APPLICATION_CORE_SOURCES} ${KRATOS_VISCOSITY_MODULATOR_TESTING_SOURCES} )
target_link_libraries( KratosViscosityModulatorCore PUBLIC KratosCore)
set_target_properties( KratosViscosityModulatorCore PROPERTIES COMPILE_DEFINITIONS "VISCOSITY_MODULATOR_APPLICATION=EXPORT,API")

## ViscosityModulatorApplication python module
pybind11_add_module( KratosViscosityModulatorApplication MODULE THIN_LTO ${KRATOS_VISCOSITY_MODULATOR_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosViscosityModulatorApplication PRIVATE KratosViscosityModulatorCore)
set_target_properties( KratosViscosityModulatorApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosViscosityModulatorCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosViscosityModulatorApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosViscosityModulatorApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosViscosityModulatorApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ViscosityModulatorApplication.py KratosMultiphysics/ViscosityModulatorApplication/__init__.py )

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
install(TARGETS KratosViscosityModulatorCore DESTINATION libs )
install(TARGETS KratosViscosityModulatorApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosViscosityModulatorCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosViscosityModulatorApplication" PARENT_SCOPE)
