set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPhaseFieldApplication ****")

include_directories(${CMAKE_SOURCE_DIR}/kratos)
include_directories(${CMAKE_SOURCE_DIR}/applications)

add_subdirectory(custom_utilities/eig)

## generate variables with the sources
set( KRATOS_PHASE_FIELD_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/phase_field_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/phase_field_fracture.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/phase_field_fracture_hybrid.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/phase_field_fracture_2.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/phase_field_kinematic_linear.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/phase_field_python_application.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPhaseFieldApplication SHARED ${KRATOS_PHASE_FIELD_APPLICATION_SOURCES})
target_link_libraries(KratosPhaseFieldApplication KratosCore)
target_link_libraries(KratosPhaseFieldApplication KratosStructuralApplication)
target_link_libraries(KratosPhaseFieldApplication eig)
target_link_libraries(KratosPhaseFieldApplication ${LAPACK_LIBRARIES})
set_target_properties(KratosPhaseFieldApplication PROPERTIES PREFIX "")
install(TARGETS KratosPhaseFieldApplication DESTINATION libs )
#exec_program(cp ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/*.so libs)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPhaseFieldApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PhaseFieldApplication.py" DESTINATION KratosMultiphysics )
