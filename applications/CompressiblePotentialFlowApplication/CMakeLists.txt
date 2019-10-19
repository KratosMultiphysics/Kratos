set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCompressiblePotentialFlowApplication ****")

################### PYBIND11
include(pybind11Tools)

if(NOT FLUID_DYNAMICS_APPLICATION )
    message(FATAL_ERROR
        "FluidDynamicsApplication is not compiled and CompressiblePotentialFlowApplication depends on it!\n"
        "  - Add '-DFLUID_DYNAMICS_APPLICATION=ON' to your build script")
endif()

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/CompressiblePotentialFlowApplication )

# generate variables with the sources
set( KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/compressible_potential_flow_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/adjoint_potential_wall_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/potential_wall_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/adjoint_analytical_incompressible_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/adjoint_base_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/adjoint_finite_difference_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/compressible_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_compressible_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/embedded_incompressible_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/incompressible_potential_flow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/apply_far_field_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_embedded_lift_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_nodal_potential_flow_velocity_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/define_2d_wake_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/define_embedded_wake_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/move_model_part_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_potential_response_function.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_response_functions/adjoint_lift_response_function_coordinates_jump.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/potential_flow_utilities.cpp
)

set( KRATOS_COMPRESSIBLE_POTENTIAL_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/compressible_potential_flow_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_response_functions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
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

if(USE_COTIRE MATCHES ON)
    cotire(KratosCompressiblePotentialFlowCore)
    cotire(KratosCompressiblePotentialFlowApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosCompressiblePotentialFlowCore DESTINATION libs )
install(TARGETS KratosCompressiblePotentialFlowApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CompressiblePotentialFlowApplication.py" DESTINATION "KratosMultiphysics/CompressiblePotentialFlowApplication" RENAME "__init__.py")

# Install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Kratos Testing.
if(${INSTALL_TESTING_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)
