set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMeshlessApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_MESHLESS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/meshless_application.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/ULF_MLSparticle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/meshless_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/meshless_application_variables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition_x.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition_y.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointforce2D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp


)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMeshlessApplication SHARED ${KRATOS_MESHLESS_APPLICATION_SOURCES})
target_link_libraries(KratosMeshlessApplication KratosCore )
set_target_properties(KratosMeshlessApplication PROPERTIES PREFIX "")
install(TARGETS KratosMeshlessApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMeshlessApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MeshlessApplication.py" DESTINATION KratosMultiphysics )
