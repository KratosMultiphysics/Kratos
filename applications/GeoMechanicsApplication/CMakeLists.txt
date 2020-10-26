set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosGeoMechanicsApplication ****")

# PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/GeoMechanicsApplication )

## Define application core sources (pure C++ parts)
set( KRATOS_GEO_MECHANICS_APPLICATION_CORE_SOURCES

    ## Application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/geo_mechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geo_mechanics_application_variables.cpp

    # Utilities
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/structural_mechanics_element_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/static_condensation_utility.cpp

    # Processes

    # Constitutive Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/elastic_isotropic_K0_3d_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_strain_K0_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_udsm_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_udsm_2D_plane_strain_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_umat_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_umat_2D_plane_strain_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_udsm_2D_interface_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_udsm_3D_interface_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_umat_2D_interface_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/small_strain_umat_3D_interface_law.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_base_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_link_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_FIC_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_strain_U_Pw_diff_order_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/drained_U_Pw_small_strain_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/undrained_U_Pw_small_strain_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_updated_lagrangian_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_U_Pw_diff_order_element.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_structural_base_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_cable_element_3D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_cr_beam_element_2D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_cr_beam_element_3D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_cr_beam_element_linear_2D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_cr_beam_element_linear_3D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_truss_element_3D2N.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/geo_truss_element_linear_3D2N.cpp

    # Conditions
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_force_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_face_load_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_FIC_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_face_load_interface_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/U_Pw_normal_flux_interface_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/general_U_Pw_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_load_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_normal_load_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_normal_fluid_flux_2D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_3D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_normal_load_3D_diff_order_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_normal_fluid_flux_3D_diff_order_condition.cpp

    # Python application

    # Strategies
)

if(${USE_PORO_MPI} MATCHES ON)

    include_directories( ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication )
    include_directories(${TRILINOS_INCLUDE_DIR})

    ## Define python interface sources
    set( KRATOS_GEO_MECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mpi_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/geo_mechanics_mpi_python_application.cpp
    )

    ## GeoMechanicsApplication core library (C++ parts)
    add_library( KratosGeoMechanicsCore SHARED ${KRATOS_GEO_MECHANICS_APPLICATION_CORE_SOURCES})
    target_link_libraries( KratosGeoMechanicsCore PUBLIC KratosCore KratosTrilinosApplication KratosStructuralMechanicsCore)


else()
    ## Define python interface sources
    set( KRATOS_GEO_MECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/geo_mechanics_python_application.cpp
    )

    ## GeoMechanicsApplication core library (C++ parts)
    add_library( KratosGeoMechanicsCore SHARED ${KRATOS_GEO_MECHANICS_APPLICATION_CORE_SOURCES})
    target_link_libraries( KratosGeoMechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)

endif()

set_target_properties( KratosGeoMechanicsCore PROPERTIES COMPILE_DEFINITIONS "GEO_MECHANICS_APPLICATION=EXPORT,API")

# define library Kratos which defines the basic python interface
## GeoMechanicsApplication python module
pybind11_add_module(KratosGeoMechanicsApplication MODULE THIN_LTO ${KRATOS_GEO_MECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosGeoMechanicsApplication PRIVATE KratosGeoMechanicsCore)
set_target_properties(KratosGeoMechanicsApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosGeoMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosGeoMechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")


if(USE_COTIRE MATCHES ON)
    cotire(KratosGeoMechanicsCore)
    cotire(KratosGeoMechanicsApplication)
endif(USE_COTIRE MATCHES ON)

## installing the resulting libraries
# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/GeoMechanicsApplication.py KratosMultiphysics/GeoMechanicsApplication/__init__.py )

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

install(TARGETS KratosGeoMechanicsCore DESTINATION libs )
install(TARGETS KratosGeoMechanicsApplication DESTINATION libs )