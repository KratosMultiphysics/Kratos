set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStructuralMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

## generate variables with the sources
set( KRATOS_STRUCTURAL_MECHANICS_APPLICATION_SOURCES
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/structural_mechanics_application_variables.cpp

  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/structural_mechanics_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_cross_sections_to_python.cpp

  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/shell_cross_section.cpp

  ##CONSTITUTIVE LAW
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/truss_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/beam_constitutive_law.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_stress.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_strain.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/elastic_isotropic_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/axisym_elastic_isotropic.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_orthotropic_2D_law.cpp

  ## ELEMENTS
  #Adding truss element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/truss_element_3D2N.cpp
   
  # Adding the beam element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cr_beam_element_3D2N.cpp

  # Adding shells and membranes elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/membrane_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/isotropic_shell_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_element_3D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thin_element_3D4N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thin_element_3D3N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_thick_element_3D3N.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/prestress_membrane_element.cpp

  # Solid elements
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_solid_element.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_small_displacement.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_total_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/axisym_updated_lagrangian.cpp
  
  # Adding the SPRISM element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/SprismElement3D6N.cpp

  # Adding the nodal concentrated element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nodal_concentrated_element.cpp

  # Adding the spring + damper element
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spring_damper_element_3D2N.cpp  

  ## CONDITIONS
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_moment_3D_condition.cpp # Beam's moment conditions
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_torque_3D_condition.cpp # Torques conditions
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/base_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_point_load_condition.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_load_condition_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_line_load_condition_2d.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_condition_3d.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosStructuralMechanicsApplication SHARED ${KRATOS_STRUCTURAL_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosStructuralMechanicsApplication KratosCore) # KratosSolidMechanicsApplication)
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosStructuralMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosStructuralMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
