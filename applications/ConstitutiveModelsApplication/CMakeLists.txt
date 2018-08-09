set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConstitutiveModelsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

# Generate variables with the sources
set( KRATOS_CONSTITUTIVE_MODELS_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_models_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_models_application_variables.cpp

	#outfitted python laws
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_python/python_outfitted_constitutive_law.cpp

	#general constitutive laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/constitutive_3D_law.cpp

	#small strain laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/small_strain_laws/small_strain_3D_law.cpp

	#large strain laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/large_strain_laws/large_strain_3D_law.cpp

	#strain rate laws
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/strain_rate_laws/strain_rate_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/strain_rate_laws/newtonian_3D_law.cpp

	#general constitutive models
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/constitutive_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/constitutive_model_data.cpp

	#elasticity models
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/linear_elastic_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/hyper_elastic_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/hypo_elastic_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/ogden_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/mooney_rivlin_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/borja_model.cpp

	#plasticity models (template, all in headers)

	#yield criteria (template, all in headers)

	#hardening rules
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/simo_linear_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/simo_exponential_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/simo_exponential_thermal_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/johnson_cook_thermal_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/baker_johnson_cook_thermal_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/exponential_damage_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/modified_exponential_damage_hardening_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/cam_clay_hardening_rule.cpp

)

## Generate variables with the sources
set( KRATOS_CONSTITUTIVE_MODELS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/constitutive_models_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

add_library(KratosConstitutiveModelsCore SHARED ${KRATOS_CONSTITUTIVE_MODELS_APPLICATION_CORE})
target_link_libraries(KratosConstitutiveModelsCore PUBLIC KratosCore)
set_target_properties(KratosConstitutiveModelsCore PROPERTIES COMPILE_DEFINITIONS "CONSTITUTIVE_MODELS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosConstitutiveModelsApplication MODULE THIN_LTO ${KRATOS_CONSTITUTIVE_MODELS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosConstitutiveModelsApplication PRIVATE KratosConstitutiveModelsCore)
set_target_properties(KratosConstitutiveModelsApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install testing files
if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosConstitutiveModelsCore)
    cotire(KratosConstitutiveModelsApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosConstitutiveModelsCore DESTINATION libs )
install(TARGETS KratosConstitutiveModelsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveModelsApplication.py" DESTINATION KratosMultiphysics )
