set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPoromechanicsApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
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

    include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application )
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
    target_link_libraries( KratosPoromechanicsCore PUBLIC KratosCore KratosTrilinosApplication KratosSolidMechanicsCore)

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
    target_link_libraries( KratosPoromechanicsCore PUBLIC KratosCore KratosSolidMechanicsCore)

endif()

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## PoromechanicsApplication core library (C++ parts)
#add_library( KratosPoromechanicsCore SHARED ${KRATOS_POROMECHANICS_APPLICATION_CORE_SOURCES})
#target_link_libraries( KratosPoromechanicsCore PRIVATE KratosCore KratosSolidMechanicsCore)
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

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosPoromechanicsCore)
    cotire(KratosPoromechanicsApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosPoromechanicsCore DESTINATION libs )
install(TARGETS KratosPoromechanicsApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PoromechanicsApplication.py" DESTINATION KratosMultiphysics )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
