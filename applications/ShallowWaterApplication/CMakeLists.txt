set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosShallowWaterApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
set( KRATOS_SHALLOW_WATER_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/shallow_water_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/shallow_water_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shallow_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/rv_swe.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cv_swe.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/swe.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/primitive_var_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conserved_var_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/euler_prim_var_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/euler_cons_var_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/nothing_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/estimate_dt_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/post_process_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shallow_water_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/replicate_model_part_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/id_renumbering_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/apply_perturbation_function_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/apply_sinusoidal_function_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/rough_porous_layer_wetting_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/negative_height_wetting_model.cpp
)

## define python interface sources
set(KRATOS_SHALLOW_WATER_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/shallow_water_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)

## define shallow water test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_SHALLOW_WATER_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## ShallowWaterApplication core library (C++ parts)
add_library( KratosShallowWaterCore SHARED ${KRATOS_SHALLOW_WATER_APPLICATION_SOURCES} ${KRATOS_SHALLOW_WATER_TESTING})
target_link_libraries( KratosShallowWaterCore PUBLIC KratosCore)
set_target_properties( KratosShallowWaterCore PROPERTIES COMPILE_DEFINITIONS "SHALLOW_WATER_APPLICATION=EXPORT,API")

## define library Kratos which defines the basic python interface
pybind11_add_module( KratosShallowWaterApplication MODULE THIN_LTO ${KRATOS_SHALLOW_WATER_APPLICATION_PYTHON_INTERFACE})
target_link_libraries( KratosShallowWaterApplication PRIVATE KratosShallowWaterCore)
set_target_properties( KratosShallowWaterApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosShallowWaterApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

install( TARGETS KratosShallowWaterCore DESTINATION libs )
install( TARGETS KratosShallowWaterApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(${INSTALL_TESTING_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ShallowWaterApplication.py" DESTINATION "KratosMultiphysics/ShallowWaterApplication" RENAME "__init__.py")
