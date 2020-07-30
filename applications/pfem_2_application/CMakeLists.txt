set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPFEM2Application ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_QUIET_WATER_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}//pfem_2_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}//pfem_2_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_2_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step_pfem_2_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/fractional_step_pfem_2_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_2fluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_2fluid_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nonewtonian_2fluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nonewtonian_2fluid_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_2fluid_2d_partintegration.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/monolithic_2fluid_3d_partintegration.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vel_enriched_2fluid_2d.cpp
	#${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/vel_enriched_2fluid_2d_nopressure.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qfluid_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qfluid_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fixed_velocity_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fixed_velocity_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fixed_pressure_2d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/fixed_pressure_3d.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/autoslip_inlet_3d.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPFEM2Application MODULE THIN_LTO ${KRATOS_QUIET_WATER_APPLICATION_SOURCES})
#add_library(KratosPFEM2Application SHARED ${KratosPFEM2Application})
target_link_libraries(KratosPFEM2Application PUBLIC KratosCore )
set_target_properties(KratosPFEM2Application PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPFEM2Application PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPFEM2Application PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosPFEM2Application)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PFEM2Application.py KratosMultiphysics/PFEM2Application/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/PFEM2Application )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/PFEM2Application
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosPFEM2Application DESTINATION libs )
