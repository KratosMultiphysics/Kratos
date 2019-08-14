set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring KratosContactStructuralMechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ContactStructuralMechanicsApplication )

## generate variables with the sources
set( KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_CORE
    ## MAIN FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application_variables.cpp

    ## UTILITIES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/derivatives_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mortar_explicit_contribution_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_preprocess.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/active_set_utilities.cpp

    ## FRICTIONAL LAWS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_frictional_laws/frictional_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_frictional_laws/frictional_law_with_derivative.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_frictional_laws/tresca_frictional_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_frictional_laws/coulomb_frictional_law.cpp

    ## PROCESSES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/aalm_adapt_penalty_value_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/master_slave_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/alm_fast_init_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/alm_variables_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_dynamic_factor_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/contact_spr_error_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/base_contact_search_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/contact_search_wrapper_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/advanced_contact_search_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/simple_contact_search_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/find_intersected_geometrical_objects_with_obb_for_contact_search_process.cpp

    ## CONDITIONS
    # Mortar condition
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/paired_condition.cpp # Paired condition
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mesh_tying_mortar_condition.cpp # Mesh tying
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mortar_contact_condition.cpp # Base class
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictionless_mortar_contact_condition.cpp # Frictionless
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictionless_components_mortar_contact_condition.cpp # Frictionless (components)
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictionless_mortar_contact_axisym_condition.cpp # Frictionless axisym
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/penalty_frictionless_mortar_contact_condition.cpp # Frictionless (penalty)
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/penalty_frictionless_mortar_contact_axisym_condition.cpp # Frictionless (penalty) axisym
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictional_mortar_contact_condition.cpp # Frictional
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictional_mortar_contact_axisym_condition.cpp # Frictional axisym
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/penalty_frictional_mortar_contact_condition.cpp # Frictional (penalty)
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/penalty_frictional_mortar_contact_axisym_condition.cpp # Frictional (penalty) axisym

    ## EXTERNAL LIBRARIES

)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## generate variables with the sources
set( KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE
    ## CUSTOM PYTHON
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_structural_mechanics_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_linear_solvers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_frictional_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/process_factory_utility.cpp
)

add_library(KratosContactStructuralMechanicsCore SHARED ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_CORE} ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosContactStructuralMechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosContactStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "CONTACT_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosContactStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosContactStructuralMechanicsApplication PRIVATE KratosContactStructuralMechanicsCore)
set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES PREFIX "")

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

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

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    ## The interface_preprocess is not compatible with cotire (solve problem with explicit instantations)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_preprocess.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    cotire(KratosContactStructuralMechanicsCore)
    cotire(KratosContactStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosContactStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosContactStructuralMechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ContactStructuralMechanicsApplication.py" DESTINATION "KratosMultiphysics/ContactStructuralMechanicsApplication" RENAME "__init__.py")
