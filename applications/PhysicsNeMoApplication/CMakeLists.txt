set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPhysicsNeMoApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## PhysicsNeMo Core sources
## NOTE: torch/physicsnemo are pure runtime Python dependencies resolved by lazy
## imports in python_scripts/ — there is deliberately no find_package(Torch), no
## CUDA detection and no USE_PHYSICSNEMO-style CMake option here.
file(GLOB_RECURSE KRATOS_PHYSICS_NEMO_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/physics_nemo_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/physics_nemo_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## PhysicsNeMo testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_PHYSICS_NEMO_APPLICATION_TESTING_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp
    )
    if(KRATOS_PHYSICS_NEMO_APPLICATION_TESTING_SOURCES)
        kratos_add_gtests(TARGET KratosPhysicsNeMoCore SOURCES "${KRATOS_PHYSICS_NEMO_APPLICATION_TESTING_SOURCES}")
    endif()
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## PhysicsNeMo python interface sources
file(GLOB_RECURSE KRATOS_PHYSICS_NEMO_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosPhysicsNeMoCore SHARED ${KRATOS_PHYSICS_NEMO_APPLICATION_CORE})
target_link_libraries(KratosPhysicsNeMoCore PUBLIC KratosCore)
set_target_properties(KratosPhysicsNeMoCore PROPERTIES COMPILE_DEFINITIONS "PHYSICS_NEMO_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPhysicsNeMoApplication MODULE THIN_LTO ${KRATOS_PHYSICS_NEMO_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPhysicsNeMoApplication PRIVATE KratosPhysicsNeMoCore)
set_target_properties(KratosPhysicsNeMoApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosPhysicsNeMoCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosPhysicsNeMoApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosPhysicsNeMoApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosPhysicsNeMoApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PhysicsNeMoApplication.py KratosMultiphysics/PhysicsNeMoApplication/__init__.py )

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
install(TARGETS KratosPhysicsNeMoCore DESTINATION libs )
install(TARGETS KratosPhysicsNeMoApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosPhysicsNeMoCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosPhysicsNeMoApplication" PARENT_SCOPE)
