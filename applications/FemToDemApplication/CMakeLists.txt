set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFemToDemApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/SolidMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/PfemFluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/MeshingApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DemStructuresCouplingApplication)

################### PYBIND11
include(pybind11Tools)

include_directories(${CMAKE_SOURCE_DIR}/kratos )
include_directories(${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication)
include_directories(${CMAKE_SOURCE_DIR}/applications/DEMApplication)
include_directories(${CMAKE_SOURCE_DIR}/applications/MeshingApplication)
include_directories(${CMAKE_SOURCE_DIR}/applications/PfemFluidDynamicsApplication)
include_directories(${CMAKE_SOURCE_DIR}/applications/DelaunayMeshingApplication)
include_directories(${CMAKE_SOURCE_DIR}/applications/DemStructuresCouplingApplication)
#include_directories(${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication)

# generate variables with the sources
set( KRATOS_FEM_TO_DEM_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/fem_to_dem_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/fem_to_dem_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fem_to_dem_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/generic_small_strain_femdem_element.cpp
    #${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/generic_large_displacement_femdem_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/extend_pressure_condition_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/update_pressure_value_pfem_conditions_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/update_pressure_volume_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/assign_pressure_id_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/expand_wet_nodes_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_initial_volume_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/regenerate_pfem_pressure_conditions_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/generate_dem_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/generate_initial_skin_DEM_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/update_dem_kinematics_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/transfer_nodal_forces_to_fem.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_normalized_free_energy_on_nodes_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/compute_sand_production.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/fix_free_velocity_on_nodes_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/remove_alone_DEM_elements_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/constitutive_law_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/FEMDEM_coupling_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/update_flag_no_remesh_femdem_boundary_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/base_solid_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/generic_total_lagrangian_femdem_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_strain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_stress.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/elastic_isotropic_3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.cpp
)

# define library Kratos which defines the basic python interface
pybind11_add_module(KratosFemToDemApplication MODULE ${KRATOS_FEM_TO_DEM_APPLICATION_SOURCES})

target_link_libraries(KratosFemToDemApplication PRIVATE KratosSolidMechanicsCore KratosDEMCore KratosMeshingCore KratosCore KratosPfemFluidDynamicsCore)
set_target_properties(KratosFemToDemApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosFemToDemApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosFemToDemApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosFemToDemApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FemToDemApplication.py KratosMultiphysics/FemToDemApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

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

# Install targets
install(TARGETS KratosFemToDemApplication DESTINATION libs )
