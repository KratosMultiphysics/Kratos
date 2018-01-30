set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosContactStructuralMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ContactStructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

## generate variables with the sources
set( KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_SOURCES
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application_variables.cpp

  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_structural_mechanics_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mappers_to_python.cpp

  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/interface_preprocess.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/tree_contact_search.cpp
  
  ## PROCESSES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/aalm_adapt_penalty_value_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/master_slave_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/alm_fast_init_process.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/alm_variables_calculation_process.cpp
  
  ## EXTERNAL LIBRARIES
  
  ## CONDITIONS
  # Mortar condition
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/paired_condition.cpp # Mesh tying
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mesh_tying_mortar_condition.cpp # Mesh tying
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_mortar_contact_condition.cpp # Base class
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictionless_mortar_contact_condition.cpp # Frictionless
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictionless_mortar_contact_axisym_condition.cpp # Frictionless axisym
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictional_mortar_contact_condition.cpp # Frictional
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/ALM_frictional_mortar_contact_axisym_condition.cpp # Frictional axisym
)

if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CONTACT_STRUCTURAL_MECHANICS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosContactStructuralMechanicsApplication SHARED ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_SOURCES} ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_TESTING_SOURCES})

target_link_libraries(KratosContactStructuralMechanicsApplication KratosCore KratosStructuralMechanicsApplication)
set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "CONTACT_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosContactStructuralMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosContactStructuralMechanicsApplication DESTINATION libs )


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

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ContactStructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
