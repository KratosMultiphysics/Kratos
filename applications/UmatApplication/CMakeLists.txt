set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosUmatApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../ConstitutiveModelsApplication )

if(NOT (${CONSTITUTIVE_MODELS_APPLICATION} MATCHES "ON"))
    message(SEND_ERROR "The ConsitutiveModels Application must be compiled before this one. Set CONSTITUTIVE_MODELS_APPLICATION=ON in the configure file.")
endif(NOT (${CONSTITUTIVE_MODELS_APPLICATION} MATCHES "ON"))

# generate variables with the sources
set( KRATOS_UMAT_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/umat_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/umat_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/umat_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/umat_3D_law.cpp

	#/${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/small_strain_umat_3D_law.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_laws/large_strain_umat_3D_law.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/mises_umat.f
	#${CMAKE_CURRENT_SOURCE_DIR}/umat_models/umat_hypo.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/hypoplastic_umat.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/umat_fortran_wrapper.f
	${CMAKE_CURRENT_SOURCE_DIR}/umat_models/xit.f
	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/small_strain_umat_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/von_mises_umat_small_strain_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/hypoplastic_umat_small_strain_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/large_strain_umat_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_models/von_mises_umat_large_strain_model.cpp

)

# define library Kratos which defines the basic python interface
add_library(KratosUmatApplication SHARED ${KRATOS_UMAT_APPLICATION_SOURCES})
target_link_libraries(KratosUmatApplication KratosCore KratosConstitutiveModelsApplication)
set_target_properties(KratosUmatApplication PROPERTIES PREFIX "")
install(TARGETS KratosUmatApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosUmatApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosUmatApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/UmatApplication.py" DESTINATION KratosMultiphysics )
