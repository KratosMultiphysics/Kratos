set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring ConstitutiveLawsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications )
include_directories( ${CMAKE_SOURCE_DIR}/applications/structural_application )

## generate variables with the sources
set( EXTERNAL_CONSTITUTIVE_LAWS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/constitutive_laws_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/umat.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/mises_umat.f
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/umat_fortran_wrapper.f
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/umat/xit.f

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosConstitutiveLawsApplication SHARED ${EXTERNAL_CONSTITUTIVE_LAWS_APPLICATION_SOURCES} )
target_link_libraries( KratosConstitutiveLawsApplication KratosCore )
set_target_properties(KratosConstitutiveLawsApplication PROPERTIES PREFIX "")
install(TARGETS KratosConstitutiveLawsApplication DESTINATION libs )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")
