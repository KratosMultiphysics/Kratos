set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosContactStructuralMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ContactStructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FSIapplication )

## generate variables with the sources
set( KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_SOURCES
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/contact_structural_mechanics_application_variables.cpp

  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/contact_structural_mechanics_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mappers_to_python.cpp

  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/tree_contact_search.cpp

  ## CONDITIONS

  # Mortar condition
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mortar_contact_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mortar_contact_2D_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mortar_contact_3D_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosContactStructuralMechanicsApplication SHARED ${KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_SOURCES})
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
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ContactStructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
