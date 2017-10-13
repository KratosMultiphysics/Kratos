set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFreezingSoilApplication ****")

include_directories( ${CMAKE_SOURCE_DIR})
include_directories( ${CMAKE_SOURCE_DIR}/kratos)
include_directories( ${CMAKE_SOURCE_DIR}/applications)
include_directories( ${CMAKE_SOURCE_DIR}/applications_nonfree)

## generate variables with the sources
set( KRATOS_FREEZING_SOIL_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/freezing_soil_application.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/freezing_soil_python_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_constitutive_laws_to_python.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_heat_flux.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_heat_convection.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_heat_radiation.cpp		
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_water_flux.cpp	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/face_load_pressure.cpp
   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/soil_2phase_rigid.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/soil_3phase.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/freezing_soil.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/unfrozen_soil.cpp   
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/solid.cpp   

	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/viscoplastic_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/clay_and_sand_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/casm.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/modified_barcelona_basic_model.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/freezing_soil_elastoplastic_model.cpp
	##${CMAKE_CURRENT_SOURCE_DIR}/constitutive_laws/freezing_soil_elastoviscoplastic_model.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosFreezingSoilApplication SHARED ${KRATOS_FREEZING_SOIL_APPLICATION_SOURCES})
target_link_libraries(KratosFreezingSoilApplication KratosCore )
set_target_properties(KratosFreezingSoilApplication PROPERTIES PREFIX "")

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFreezingSoilApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFreezingSoilApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosFreezingSoilApplication DESTINATION libs )


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)



# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/FreezingSoilApplication.py" DESTINATION KratosMultiphysics )
