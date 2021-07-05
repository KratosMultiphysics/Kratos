set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConstitutiveLawsApplication ****")

################### PYBIND11
include(pybind11Tools)

# temp dependency until the remaining basic laws are moved to the Core
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## ConstitutiveLaws Core sources
file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## ConstitutiveLaws testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ConstitutiveLaws python interface sources
file(GLOB_RECURSE KRATOS_CONSTITUTIVE_LAWS_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosConstitutiveLawsCore SHARED ${KRATOS_CONSTITUTIVE_LAWS_APPLICATION_CORE} ${KRATOS_CONSTITUTIVE_LAWS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosConstitutiveLawsCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosConstitutiveLawsCore PROPERTIES COMPILE_DEFINITIONS "CONSTITUTIVE_LAWS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosConstitutiveLawsApplication MODULE THIN_LTO ${KRATOS_CONSTITUTIVE_LAWS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosConstitutiveLawsApplication PRIVATE KratosConstitutiveLawsCore)
set_target_properties(KratosConstitutiveLawsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosConstitutiveLawsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosConstitutiveLawsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosConstitutiveLawsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosConstitutiveLawsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveLawsApplication.py KratosMultiphysics/ConstitutiveLawsApplication/__init__.py )

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
install(TARGETS KratosConstitutiveLawsCore DESTINATION libs )
install(TARGETS KratosConstitutiveLawsApplication DESTINATION libs )
