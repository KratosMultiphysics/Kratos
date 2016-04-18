set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPoromechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )

## generate variables with the sources
set( KRATOS_POROMECHANICS_APPLICATION_SOURCES

	${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/poromechanics_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_2D_law.cpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_element.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_interface_element.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_link_interface_element.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_FIC_element.cpp #TODO
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_strain_U_Pw_diff_order_element.cpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_force_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_FIC_condition.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_face_load_interface_condition.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_interface_condition.cpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/general_U_Pw_diff_order_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_load_2D_diff_order_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_normal_load_2D_diff_order_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_3D_diff_order_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_normal_load_3D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPoromechanicsApplication SHARED ${KRATOS_POROMECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosPoromechanicsApplication KratosCore KratosSolidMechanicsApplication)
set_target_properties(KratosPoromechanicsApplication PROPERTIES PREFIX "")
set_target_properties(KratosPoromechanicsApplication PROPERTIES COMPILE_DEFINITIONS "POROMECHANICS_APPLICATION=EXPORT,API")
install(TARGETS KratosPoromechanicsApplication DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PoromechanicsApplication.py" DESTINATION KratosMultiphysics )
