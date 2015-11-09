set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPFEMFluidDynamicApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_PFEM_FLUID_DYNAMIC_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//pfem_fluid_dynamic_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}//pfem_fluid_dynamic_application_variables.cpp 

	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_fluid_dynamic_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPFEMFluidDynamicApplication SHARED ${KRATOS_PFEM_FLUID_DYNAMIC_APPLICATION_SOURCES})
target_link_libraries(KratosPFEMFluidDynamicApplication KratosCore )
set_target_properties(KratosPFEMFluidDynamicApplication PROPERTIES PREFIX "")
install(TARGETS KratosPFEMFluidDynamicApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPFEMFluidDynamicApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPFEMFluidDynamicApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PFEMFluidDynamicApplication.py" DESTINATION KratosMultiphysics )
