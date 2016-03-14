set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDEMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_DEM_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/DEM_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/DEM_python_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp 
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spheric_particle.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/nanoparticle.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/spheric_continuum_particle.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/thermal_spheric_particle.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/Particle_Contact_Element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cylinder_particle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cylinder_continuum_particle.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/linecluster3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/RigidFace.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/RigidEdge.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dem_wall.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/mapping_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_discontinuum_constitutive_law.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_continuum_constitutive_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/cubecluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/pillcluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ellipsoidcluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ringcluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/cuboidcluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/cornkernelcluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/corn3cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/soybeancluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/rock1cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/rock2cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/wheat5cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast1cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast2cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast3cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast4cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/capsulecluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast5cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/custom_clusters/ballast6cluster3D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/forward_euler_scheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/symplectic_euler_scheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/newmark_beta_scheme.cpp        
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/dem_integration_scheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/verlet_velocity_scheme.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/schemes/symplectic_euler_scheme.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_2D_CL.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_Dempack_torque_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_KDEM_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_ExponentialHC_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_viscous_Coulomb_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_viscous_Coulomb_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Bentonite_Colloid_CL.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_JKR_cohesive_law.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_DMT_cohesive_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Linear_viscous_Coulomb_2D_CL.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/DEM_D_Hertz_viscous_Coulomb_2D_CL.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/inlet.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/create_and_destroy.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/properties_proxies.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/dem_fem_utilities.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosDEMApplication SHARED ${KRATOS_DEM_APPLICATION_SOURCES})
#if(NOT EXISTS "${PROJECT_SOURCE_DIR}/applications/DEM_application/kratos_DEMApplication_export_dll.h")
#GENERATE_EXPORT_HEADER( KratosDEMApplication
#             BASE_NAME KratosDEMApplication
#             EXPORT_MACRO_NAME KRATOS_DEMAPPLICATION_EXPORT_DLL
#             EXPORT_FILE_NAME ../../../applications/DEM_application/kratos_DEMApplication_export_dll.h
#             STATIC_DEFINE KratosDemApplication_BUILT_AS_STATIC
#)
#endif(NOT EXISTS "${PROJECT_SOURCE_DIR}/applications/DEM_application/kratos_DEMApplication_export_dll.h")
ADD_COMPILER_EXPORT_FLAGS()
target_link_libraries(KratosDEMApplication KratosCore )
set_target_properties(KratosDEMApplication PROPERTIES PREFIX "")
install(TARGETS KratosDEMApplication DESTINATION libs )

if(${ACTIVATE_DEBUG_MACRO} MATCHES ON)        #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDEMApplication PROPERTIES SUFFIX .pyd)
    set_target_properties(KratosDEMApplication PROPERTIES COMPILE_DEFINITIONS "DEM_APPLICATION=EXPORT,API")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDEMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosDEMApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/DEMApplication.py" DESTINATION KratosMultiphysics )
