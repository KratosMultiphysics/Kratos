set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIGAStructuralMechanicsApplication ****")

include(pybind11Tools)

if( NOT DEFINED STRUCTURAL_MECHANICS_APPLICATION )
	message( FATAL_ERROR "STRUCTURAL_MECHANICS_APPLICATION not defined" )
endif()

if( STRUCTURAL_MECHANICS_APPLICATION MATCHES OFF )
	message( FATAL_ERROR "STRUCTURAL_MECHANICS_APPLICATION not being compiled" )
endif()

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
#include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

# generate variables with the sources
set( KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_SOURCES
	## MAIN FILES
	${CMAKE_CURRENT_SOURCE_DIR}/iga_structural_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/iga_structural_mechanics_application_variables.cpp


	## UTILITIES
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/NurbsShapeFunctionModeler.cpp

	##ELEMENTS
	# Adding the discrete base elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_discrete_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/surface_base_discrete_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/curve_base_discrete_element.cpp

	# Adding the discrete elements
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_kl_discrete_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/membrane_discrete_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/truss_discrete_element.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_base_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_base_surface_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_shell_kl_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_shell_kl_thick_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_membrane_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_shell_element.cpp

	##CONDITIONS
	# Adding the discrete base conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/base_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_base_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/curve_base_discrete_condition.cpp

	# Adding the discrete conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_base_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_surface_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_curve_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/load_point_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/support_penalty_curve_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/support_strong_discrete_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/support_penalty_point_discrete_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_base_condition.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_force_interface_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_support_rotation_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_surface_support_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_base_coupling_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_lagrange_coupling_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_lagrange_coupling_condition_2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_penalty_coupling_rotation_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_penalty_coupling_crack_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/iga_flags.cpp


	# ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/damage_tc_plane_stress_2d_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/plane_stress_tc_damage_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/plane_stress_2d_tc_damage_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/plane_stress_2d_kinematically_enriched_law.cpp
)
## generate variables with the sources
set( KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE
	## CUSTOM PYTHON
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/iga_structural_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)
# define library Kratos which defines the basic python interface
#add_library(KratosIGAStructuralMechanicsCore SHARED ${KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_SOURCES})
#target_link_libraries( KratosIGAStructuralMechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)
#set_target_properties(KratosIGAStructuralMechanicsCore PROPERTIES COMPILE_DEFINITIONS "IGA_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")

pybind11_add_module(KratosIGAStructuralMechanicsApplication MODULE THIN_LTO ${KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_SOURCES}
                                                                            ${KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosIGAStructuralMechanicsApplication PUBLIC KratosCore KratosStructuralMechanicsCore )
set_target_properties(KratosIGAStructuralMechanicsApplication PROPERTIES COMPILE_DEFINITIONS "IGA_STRUCTURAL_MECHANICS_APPLICATION=EXPORT,API")
set_target_properties(KratosIGAStructuralMechanicsApplication PROPERTIES PREFIX "")



# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosIGAStructuralMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosIGAStructuralMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
	get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    cotire(KratosIGAStructuralMechanicsCore)
    cotire(KratosIGAStructuralMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

#install(TARGETS KratosIGAStructuralMechanicsCore DESTINATION libs )
install(TARGETS KratosIGAStructuralMechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/IGAStructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
