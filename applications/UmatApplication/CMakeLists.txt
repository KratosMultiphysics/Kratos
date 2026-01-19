set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosUmatApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ConstitutiveModelsApplication )

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ConstitutiveModelsApplication)

# Enable Fortran compilation
enable_language(Fortran)

# Generate variables with the sources
set( KRATOS_UMAT_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/umat_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/umat_application_variables.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/mises_umat.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/hypoplastic_umat.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/NiceNames.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/TensorModule.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/CompLib_new.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/umat_fortran_wrapper.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/xit.f

	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/small_strain_umat_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/von_mises_umat_small_strain_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/fabric_umat_small_strain_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/hypoplastic_umat_small_strain_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/large_strain_umat_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/von_mises_umat_large_strain_model.cpp
)

## Generate variables with the sources
set( KRATOS_UMAT_APPLICATION_PYTHON_INTERFACE
    # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/umat_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

add_library(KratosUmatCore SHARED ${KRATOS_UMAT_APPLICATION_CORE})
target_link_libraries(KratosUmatCore PUBLIC KratosCore KratosConstitutiveModelsCore)
set_target_properties(KratosUmatCore PROPERTIES COMPILE_DEFINITIONS "UMAT_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosUmatApplication MODULE THIN_LTO ${KRATOS_UMAT_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosUmatApplication PRIVATE KratosUmatCore)
set_target_properties(KratosUmatApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosUmatApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosUmatApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosUmatApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/UmatApplication.py KratosMultiphysics/UmatApplication/__init__.py )

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
install(TARGETS KratosUmatCore DESTINATION libs )
install(TARGETS KratosUmatApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosUmatCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosUmatApplication" PARENT_SCOPE)