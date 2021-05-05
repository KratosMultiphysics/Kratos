set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPlasmaDynamicsApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/SwimmingDEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ConvectionDiffusionApplication)

include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Dependecies on other applications
include_directories( ${KRATOS_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/ConvectionDiffusionApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/SwimmingDEMApplication )
# include_directories( ${KRATOS_SOURCE_DIR}/applications/PlasmaDynamicsApplication )

## Generate variables with the sources
set( KRATOS_PLASMA_DYNAMICS_APPLICATION_CORE_SOURCES
    ## Application class
    ${CMAKE_CURRENT_SOURCE_DIR}/plasma_dynamics_application.cpp


)

## Define python interface sources
set( KRATOS_PLASMA_DYNAMICS_APPLICATION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/plasma_dynamics_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPlasmaDynamicsCore SHARED ${KRATOS_PLASMA_DYNAMICS_APPLICATION_CORE_SOURCES})
target_link_libraries(KratosPlasmaDynamicsCore KratosDEMCore KratosFluidDynamicsCore KratosSwimmingDEMCore KratosConvectionDiffusionCore KratosCore)
set_target_properties(KratosPlasmaDynamicsCore PROPERTIES COMPILE_DEFINITIONS "PLASMA_DYNAMICS_APPLICATION=EXPORT,API")

## python module
pybind11_add_module(KratosPlasmaDynamicsApplication MODULE ${KRATOS_PLASMA_DYNAMICS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosPlasmaDynamicsApplication PUBLIC KratosPlasmaDynamicsCore)
set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPlasmaDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosPlasmaDynamicsCore)
  cotire(KratosPlasmaDynamicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PlasmaDynamicsApplication.py KratosMultiphysics/PlasmaDynamicsApplication/__init__.py )

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

## Install targets
install(TARGETS KratosPlasmaDynamicsCore DESTINATION libs )
install(TARGETS KratosPlasmaDynamicsApplication DESTINATION libs )
