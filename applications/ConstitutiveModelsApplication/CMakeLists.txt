set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConstitutiveModelsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/elasticity_models/tamagnini_model.cpp

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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/plasticity_models/hardening_rules/gens_nova_hardening_rule.cpp

    #utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/properties_layout.cpp

	#processes
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/non_local_plasticity_process.cpp

)

## Generate variables with the sources
set( KRATOS_CONSTITUTIVE_MODELS_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/constitutive_models_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

add_library(KratosConstitutiveModelsCore SHARED ${KRATOS_CONSTITUTIVE_MODELS_APPLICATION_CORE})
target_link_libraries(KratosConstitutiveModelsCore PUBLIC KratosCore)
set_target_properties(KratosConstitutiveModelsCore PROPERTIES COMPILE_DEFINITIONS "CONSTITUTIVE_MODELS_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosConstitutiveModelsApplication MODULE THIN_LTO ${KRATOS_CONSTITUTIVE_MODELS_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosConstitutiveModelsApplication PRIVATE KratosConstitutiveModelsCore)
set_target_properties(KratosConstitutiveModelsApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosConstitutiveModelsCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosConstitutiveModelsApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosConstitutiveModelsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ConstitutiveModelsApplication.py KratosMultiphysics/ConstitutiveModelsApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosConstitutiveModelsCore DESTINATION libs )
install(TARGETS KratosConstitutiveModelsApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosConstitutiveModelsCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosConstitutiveModelsApplication" PARENT_SCOPE)