set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosChimeraApplication****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

## ChimeraApplication Core sources
file(GLOB_RECURSE KRATOS_CHIMERA_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/chimera_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimera_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## ChimeraApplication testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CHIMERA_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## ChimeraApplication python interface sources
file(GLOB_RECURSE KRATOS_CHIMERA_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosChimeraApplicationCore SHARED ${KRATOS_CHIMERA_APPLICATION_CORE} ${KRATOS_CHIMERA_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosChimeraApplicationCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties(KratosChimeraApplicationCore PROPERTIES COMPILE_DEFINITIONS "CHIMERA_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosChimeraApplication MODULE THIN_LTO ${KRATOS_CHIMERA_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosChimeraApplication PRIVATE KratosChimeraApplicationCore)
set_target_properties(KratosChimeraApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosChimeraApplicationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosChimeraApplicationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosChimeraApplicationCore)
    cotire(KratosChimeraApplicationApplication)
endif(USE_COTIRE MATCHES ON)



# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ChimeraApplication.py KratosMultiphysics/ChimeraApplication/__init__.py )

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

install(TARGETS KratosChimeraApplicationCore DESTINATION libs )
install(TARGETS KratosChimeraApplication DESTINATION libs )
