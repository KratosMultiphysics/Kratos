set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosThermalDEMApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication)

## Dependecies on other applications
include(pybind11Tools)
include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/ThermalDEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/DelaunayMeshingApplication )

## Set a default value for triangle and tetgen in case the user does not set it
## to avoid problems with the define linkage block. By default we will compile the lib
if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_DIR "${KRATOS_SOURCE_DIR}/external_libraries/triangle")
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../DelaunayMeshingApplication/external_libraries/tetgen")
endif(NOT DEFINED ${INCLUDE_TETGEN})

set(MESHER_LIBRARIES "")
if(${INCLUDE_TRIANGLE} MATCHES ON)
  include_directories(${TRIANGLE_DIR})
  list(APPEND MESHER_LIBRARIES triangle)
endif(${INCLUDE_TRIANGLE} MATCHES ON)

if(${INCLUDE_TETGEN} MATCHES ON)
  include_directories(${TETGEN_DIR})
  list(APPEND MESHER_LIBRARIES tetgen)
endif(${INCLUDE_TETGEN} MATCHES ON)

if(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  list(APPEND MESHER_LIBRARIES m)
endif(NOT ${CMAKE_SYSTEM_NAME} MATCHES "Windows")

## Generate variables with the sources
set( KRATOS_THERMAL_DEM_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/thermal_dem_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/thermal_dem_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/heat_exchange_mechanism.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/heat_generation_mechanism.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/sintering_continuum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_bob_complete.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_bob_modified.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_bob_simple.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_collision.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/direct_conduction_pipe.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/indirect_conduction_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/indirect_conduction_surround_layer.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/indirect_conduction_vargas.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/indirect_conduction_voronoi_a.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/conduction/indirect_conduction_voronoi_b.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/convection/convection_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/convection/nusselt_gunn.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/convection/nusselt_hanz_marshall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/convection/nusselt_li_mason.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/convection/nusselt_whitaker.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generation/generation_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/generation/generation_dissipation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/radiation/radiation_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/radiation/radiation_continuum_krause.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/radiation/radiation_continuum_zhou.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_model.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_lu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_zhou.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_morris_area.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_morris_area_time.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_rangel_area.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/real_contact/real_contact_rangel_area_time.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_spheric_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_spheric_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/sintering_spheric_continuum_particle.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/strategies/thermal_explicit_solver_strategy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/thermal_dem_integration_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/thermal_forward_euler_scheme.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/set_thermal_data_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/numerical_integration_method.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/numerical_integration_adaptive_simpson.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/tesselation_utilities_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/tesselation_utilities_3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/graph_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/heat_map_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rve_wall_boundary_thermal_2d.cpp
)

## Define python interface sources
set( KRATOS_THERMAL_DEM_APPLICATION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/thermal_dem_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

## Define library Kratos which defines the basic python interface
add_library(KratosThermalDEMCore SHARED ${KRATOS_THERMAL_DEM_APPLICATION_CORE_SOURCES})
target_link_libraries(KratosThermalDEMCore KratosDEMCore KratosCore ${MESHER_LIBRARIES})
set_target_properties(KratosThermalDEMCore PROPERTIES COMPILE_DEFINITIONS "THERMAL_DEM_APPLICATION=EXPORT,API")

## Python module
pybind11_add_module(KratosThermalDEMApplication MODULE ${KRATOS_THERMAL_DEM_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosThermalDEMApplication PUBLIC KratosThermalDEMCore)
set_target_properties(KratosThermalDEMApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosThermalDEMCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosThermalDEMApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Change the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosThermalDEMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# Change the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosThermalDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphysics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/ThermalDEMApplication.py KratosMultiphysics/ThermalDEMApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/tests KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

## Install targets
install(TARGETS KratosThermalDEMCore DESTINATION libs )
install(TARGETS KratosThermalDEMApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosThermalDEMCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosThermalDEMApplication" PARENT_SCOPE)
