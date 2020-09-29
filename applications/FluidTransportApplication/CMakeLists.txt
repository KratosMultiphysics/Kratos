set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidTransportApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidTransportApplication )


## Define application core sources (pure C++ parts)
set( KRATOS_FLUID_TRANSPORT_APPLICATION_CORE_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/fluid_transport_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/fluid_transport_application_variables.cpp

  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/steady_convection_diffusion_FIC_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/transient_convection_diffusion_FIC_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/transient_convection_diffusion_FIC_explicit_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/transient_convection_diffusion_PFEM2_FIC_element.cpp
)

## Define python interface sources
set( KRATOS_FLUID_TRANSPORT_APPLICATION_PYTHON_INTERFACE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fluid_transport_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)

###############################################################################
## FluidTransportApplication core library (C++ parts)
add_library(KratosFluidTransportCore SHARED ${KRATOS_FLUID_TRANSPORT_APPLICATION_CORE_SOURCES})
target_link_libraries(KratosFluidTransportCore PRIVATE KratosCore KratosFluidDynamicsCore)

## FluidTransportApplication python module
pybind11_add_module(KratosFluidTransportApplication MODULE ${KRATOS_FLUID_TRANSPORT_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosFluidTransportApplication PRIVATE KratosFluidTransportCore)
set_target_properties(KratosFluidTransportApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosFluidTransportCore)
    cotire(KratosFluidTransportApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidTransportApplication.py KratosMultiphysics/FluidTransportApplication/__init__.py )

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
install(TARGETS KratosFluidTransportCore DESTINATION libs )
install(TARGETS KratosFluidTransportApplication DESTINATION libs )