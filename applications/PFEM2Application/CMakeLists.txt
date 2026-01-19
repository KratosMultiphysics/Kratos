set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPFEM2Application ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

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

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosPFEM2Application PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPFEM2Application PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPFEM2Application PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PFEM2Application.py KratosMultiphysics/PFEM2Application/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/PFEM2Application )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
	kratos_python_install_tests(
        DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests"
        DESTINATION "applications/${CURRENT_DIR_NAME}" 
        INSTALL_PYTHON_USING_LINKS "${INSTALL_PYTHON_USING_LINKS}"
        EXCLUSION_PATTERNS "\.git$" "\.c$" "\.h$" "\.cpp$" "\.hpp$"
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosPFEM2Application DESTINATION libs )

# Define custom targets
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosPFEM2Application" PARENT_SCOPE)