set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPoromechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/PoromechanicsApplication )

## Define application core sources (pure C++ parts)
set( KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES

    ## Application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/poromechanics_application_variables.cpp

    # Utilities

    # Processes

    # Constitutive Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/local_damage_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/nonlocal_damage_flow_rule.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/bilinear_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/exponential_cohesive_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/exponential_cohesive_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/simo_ju_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/modified_mises_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/modified_mises_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/modified_mises_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/history_linear_elastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/history_linear_elastic_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/history_linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_hardening_laws/modified_exponential_damage_hardening_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_yield_criteria/modified_mises_yield_criterion.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/custom_flow_rules/isotropic_damage_flow_rule.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plane_stress_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/hyperelastic_plastic_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_plastic_3D_law.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_link_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/U_Pw_small_strain_FIC_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_strain_U_Pw_diff_order_element.cpp

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
)

if(${USE_PORO_MPI} MATCHES ON)

    include_directories( ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication )
    include_directories(${TRILINOS_INCLUDE_DIR})

    ## Define python interface sources
    set( KRATOS_POROMECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mpi_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/poromechanics_mpi_python_application.cpp
    )

    ## PoromechanicsApplication core library (C++ parts)
    add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
    target_link_libraries( KratosPoromechanicsCore PUBLIC KratosCore KratosTrilinosApplication KratosStructuralMechanicsCore)

else()

    ## Define python interface sources
    set( KRATOS_POROMECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/poromechanics_python_application.cpp
    )

    ## PoromechanicsApplication core library (C++ parts)
    add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
    target_link_libraries( KratosPoromechanicsCore PUBLIC KratosCore KratosStructuralMechanicsCore)

endif()

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## PoromechanicsApplication core library (C++ parts)
#add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
#target_link_libraries( KratosPoromechanicsCore PRIVATE KratosCore KratosStructuralMechanicsCore)
set_target_properties( KratosPoromechanicsCore PROPERTIES COMPILE_DEFINITIONS "POROMECHANICS_APPLICATION=EXPORT,API")

## PoromechanicsApplication python module
pybind11_add_module(KratosPoromechanicsApplication MODULE THIN_LTO ${KRATOS_POROMECHANICS_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosPoromechanicsApplication PRIVATE KratosPoromechanicsCore)
set_target_properties(KratosPoromechanicsApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPoromechanicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosPoromechanicsCore)
    cotire(KratosPoromechanicsApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PoromechanicsApplication.py KratosMultiphysics/PoromechanicsApplication/__init__.py )

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
install(TARGETS KratosPoromechanicsCore DESTINATION libs )
install(TARGETS KratosPoromechanicsApplication DESTINATION libs )
