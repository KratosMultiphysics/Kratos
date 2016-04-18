set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDamApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/PoromechanicsApplication )

## generate variables with the sources
set( KRATOS_DAM_APPLICATION_SOURCES   
	${CMAKE_CURRENT_SOURCE_DIR}/dam_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dam_application_variables.cpp
    
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dam_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_thermo_mechanic_element.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_interface_element.cpp
    
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_strain.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_stress.cpp
    
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosDamApplication SHARED ${KRATOS_DAM_APPLICATION_SOURCES})
target_link_libraries(KratosDamApplication KratosCore KratosSolidMechanicsApplication KratosPoromechanicsApplication)
set_target_properties(KratosDamApplication PROPERTIES PREFIX "")
install(TARGETS KratosDamApplication DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDamApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDamApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DamApplication.py" DESTINATION KratosMultiphysics )
