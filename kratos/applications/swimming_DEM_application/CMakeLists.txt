set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosSwimmingDEMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEM_application )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

## generate variables with the sources
set( KRATOS_SWIMMING_DEM_APPLICATION_SOURCES    

    ${CMAKE_CURRENT_SOURCE_DIR}/swimming_DEM_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_dem_coupled.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_dem_coupled_weak.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/monolithic_dem_coupled_wall_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/shell_rigid.cpp    
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spheric_swimming_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/swimming_DEM_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/terminal_velocity_scheme.cpp    
)


###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosSwimmingDEMApplication SHARED ${KRATOS_SWIMMING_DEM_APPLICATION_SOURCES})
target_link_libraries(KratosSwimmingDEMApplication KratosDEMApplication KratosFluidDynamicsApplication KratosCore)
set_target_properties(KratosSwimmingDEMApplication PROPERTIES PREFIX "")
install(TARGETS KratosSwimmingDEMApplication DESTINATION libs )

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .pyd)
    set_target_properties(KratosSwimmingDEMApplication PROPERTIES COMPILE_DEFINITIONS "SWIMMING_DEM_APPLICATION=EXPORT,API")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosSwimmingDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosSwimmingDEMApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/SwimmingDEMApplication.py" DESTINATION KratosMultiphysics )
