set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosIGAStructuralMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication ) 
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication ) 
include_directories( ${CMAKE_SOURCE_DIR}/applications/MeshingApplication )

# generate variables with the sources
set( KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_SOURCES
	## MAIN FILES
	${CMAKE_CURRENT_SOURCE_DIR}/iga_structural_mechanics_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/iga_structural_mechanics_application_variables.cpp
	
	## CUSTOM PYTHON
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/iga_structural_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	
	## UTILITIES
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/NurbsShapeFunctionModeler.cpp		
		
	##ELEMENTS
	# Adding the meshless base element
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_base_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_membrane_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_laplace_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_shell_element.cpp
		
	##CONDITIONS  ContinuityConditionPenalty
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_base_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_load_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_support_rotation_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_base_coupling_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_lagrange_coupling_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/meshless_penalty_coupling_rotation_condition.cpp
)

# define library Kratos which defines the basic python interface
add_library(KratosIGAStructuralMechanicsApplication SHARED ${KRATOS_IGA_STRUCTURAL_MECHANICS_APPLICATION_SOURCES})
target_link_libraries(KratosIGAStructuralMechanicsApplication KratosStructuralMechanicsApplication KratosCore )
set_target_properties(KratosIGAStructuralMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosIGAStructuralMechanicsApplication DESTINATION libs )

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

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/IGAStructuralMechanicsApplication.py" DESTINATION KratosMultiphysics )
