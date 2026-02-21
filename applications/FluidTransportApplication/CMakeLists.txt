set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFluidTransportApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/FluidTransportApplication )


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
target_link_libraries(KratosFluidTransportCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties( KratosFluidTransportCore PROPERTIES COMPILE_DEFINITIONS "FLUID_TRANSPORT_APPLICATION=EXPORT,API")

## FluidTransportApplication python module
pybind11_add_module(KratosFluidTransportApplication MODULE THIN_LTO ${KRATOS_FLUID_TRANSPORT_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosFluidTransportApplication PRIVATE KratosFluidTransportCore)
set_target_properties(KratosFluidTransportApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFluidTransportCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFluidTransportApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFluidTransportApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FluidTransportApplication.py KratosMultiphysics/FluidTransportApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosFluidTransportCore DESTINATION libs )
install(TARGETS KratosFluidTransportApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFluidTransportCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFluidTransportApplication" PARENT_SCOPE)