set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCompressiblePotentialFlowApplication ****")

################### PYBIND11
include(pybind11Tools)

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/CompressiblePotentialFlowApplication )

# generate variables with the sources
file( GLOB_RECURSE KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

file( GLOB_RECURSE KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)
## define compressible potential application test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosCompressiblePotentialFlowCore SHARED ${KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_CORE} ${KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosCompressiblePotentialFlowCore PUBLIC KratosCore KratosFluidDynamicsCore)
set_target_properties(KratosCompressiblePotentialFlowCore PROPERTIES COMPILE_DEFINITIONS "COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION=EXPORT,API")

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosCompressiblePotentialFlowApplication MODULE THIN_LTO ${KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosCompressiblePotentialFlowApplication PRIVATE KratosCompressiblePotentialFlowCore )
set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosCompressiblePotentialFlowApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosCompressiblePotentialFlowCore)
    cotire(KratosCompressiblePotentialFlowApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/CompressiblePotentialFlowApplication.py KratosMultiphysics/CompressiblePotentialFlowApplication/__init__.py )

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
install(TARGETS KratosCompressiblePotentialFlowCore DESTINATION libs )
install(TARGETS KratosCompressiblePotentialFlowApplication DESTINATION libs )
