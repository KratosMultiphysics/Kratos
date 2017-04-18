set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConstitutiveModelsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_CONSTITUTIVE_MODELS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_models_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_models_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/constitutive_models_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	
	#outfitted python laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/python_outfitted_constitutive_law.cpp

	#general constitutive laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/elastic_3D_law.cpp

	#isotropic linear elastic laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_laws/linear_elastic_3D_law.cpp

	#isotropic hyperelastic laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/hyperelastic_laws/hyperelastic_3D_law.cpp
		
	#isotropic linear elastic plasticity laws
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_3D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_strain_2D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_plane_stress_2D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_plastic_laws/linear_elastic_plastic_axysimmetric_2D_law.cpp

	#isotropic hyperelastic plasticity laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/hyperelastic_plastic_laws/hyperelastic_plastic_3D_law.cpp


	#isotropic linear elastic damage laws
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_3D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_strain_2D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_simo_ju_plane_stress_2D_law.cpp

	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_3D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_strain_2D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/linear_elastic_damage_laws/isotropic_damage_modified_mises_plane_stress_2D_law.cpp


	#general constitutive models
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/constitutive_model_data.cpp

	#elasticity models
		
	#hyperelastic models
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/hyperelastic_models/hyperelastic_model.cpp

	#plastic models (template, all in headers)
	
	#yield criteria (template, all in headers)

	#hardening laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_laws/non_linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_laws/linear_isotropic_kinematic_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_laws/exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_laws/modified_exponential_damage_hardening_law.cpp

	
)

# define library Kratos which defines the basic python interface
add_library(KratosConstitutiveModelsApplication SHARED ${KRATOS_CONSTITUTIVE_MODELS_APPLICATION_SOURCES})
target_link_libraries(KratosConstitutiveModelsApplication KratosCore )
set_target_properties(KratosConstitutiveModelsApplication PROPERTIES PREFIX "")
set_target_properties(KratosConstitutiveModelsApplication PROPERTIES COMPILE_DEFINITIONS "CONSTITUTIVE_MODELS_APPLICATION=EXPORT,API")
install(TARGETS KratosConstitutiveModelsApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveModelsApplication.py" DESTINATION KratosMultiphysics )
