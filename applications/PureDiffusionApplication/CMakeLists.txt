set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPureDiffusionApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/convection_diffusion_application )

## generate variables with the sources
set( KRATOS_PURE_DIFFUSION_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/pure_diffusion_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pure_diffusion_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pure_diffusion_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/poisson_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/projected_swe.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/point_source.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPureDiffusionApplication SHARED ${KRATOS_PURE_DIFFUSION_APPLICATION_SOURCES})
target_link_libraries(KratosPureDiffusionApplication KratosCore )
set_target_properties(KratosPureDiffusionApplication PROPERTIES PREFIX "")
install(TARGETS KratosPureDiffusionApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPureDiffusionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PureDiffusionApplication.py" DESTINATION KratosMultiphysics )
