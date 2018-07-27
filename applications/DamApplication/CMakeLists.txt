set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDamApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/PoromechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/DamApplication )

## Define application core sources (pure C++ parts)
set( KRATOS_DAM_APPLICATION_CORE_SOURCES

    ## Application class and variables
    ${CMAKE_CURRENT_SOURCE_DIR}/dam_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dam_application_variables.cpp

    # Utilities

    # Pocesses

    # Constitutive Laws
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_strain.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_stress.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_3D_law_nodal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_2D_plane_strain_nodal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_elastic_2D_plane_stress_nodal.cpp 
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_3D_law_nodal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_strain_nodal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_linear_elastic_2D_plane_stress_nodal.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_local_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_local_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_local_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_simo_ju_nonlocal_damage_plane_stress_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_modified_mises_nonlocal_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_strain_2D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/thermal_modified_mises_nonlocal_damage_plane_stress_2D_law.cpp

    # Elements
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_thermo_mechanic_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_interface_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/wave_equation_element.cpp

    # Conditions
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/UP_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/free_surface_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/infinite_domain_condition.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/added_mass_condition.cpp
)

if(${USE_DAM_MPI} MATCHES ON)

    include_directories( ${CMAKE_SOURCE_DIR}/applications/trilinos_application )
    include_directories(${TRILINOS_INCLUDE_DIR})

    ## Define python interface sources
    set( KRATOS_DAM_APPLICATION_PYTHON_INTERFACE_SOURCES

        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_mpi_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dam_mpi_python_application.cpp
    )

    ## DamApplication core library (C++ parts)
    add_library(KratosDamCore SHARED ${KRATOS_DAM_APPLICATION_CORE_SOURCES})
    target_link_libraries(KratosDamCore PUBLIC KratosCore KratosTrilinosApplication KratosSolidMechanicsCore KratosPoromechanicsCore)

else()

    ## Define python interface sources
    set( KRATOS_DAM_APPLICATION_PYTHON_INTERFACE_SOURCES
        
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dam_python_application.cpp
    )
    
    ## DamApplication core library (C++ parts)
    add_library(KratosDamCore SHARED ${KRATOS_DAM_APPLICATION_CORE_SOURCES})
    target_link_libraries(KratosDamCore PUBLIC KratosCore KratosSolidMechanicsCore KratosPoromechanicsCore)

endif()

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

###############################################################################
## DamApplication core library (C++ parts)
#add_library(KratosDamCore SHARED ${KRATOS_DAM_APPLICATION_CORE_SOURCES})
#target_link_libraries(KratosDamCore PRIVATE KratosCore KratosSolidMechanicsCore KratosPoromechanicsCore)
set_target_properties(KratosDamCore PROPERTIES COMPILE_DEFINITIONS "DAM_APPLICATION=EXPORT,API")

## DamApplication python module
pybind11_add_module(KratosDamApplication MODULE THIN_LTO ${KRATOS_DAM_APPLICATION_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosDamApplication PRIVATE KratosDamCore)
set_target_properties(KratosDamApplication PROPERTIES PREFIX "")

###############################################################################
# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosDamApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDamApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosDamCore)
    cotire(KratosDamApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosDamCore DESTINATION libs )
install(TARGETS KratosDamApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DamApplication.py" DESTINATION KratosMultiphysics )

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
