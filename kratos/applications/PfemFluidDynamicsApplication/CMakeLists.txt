set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemFluidDynamicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemFluidDynamicsApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication/custom_external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../MeshingApplication/custom_external_libraries/tetgen/tetgen1.4.3 )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)


## generate variables with the sources
set( KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_solid_mechanics_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_fluid_dynamics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemFluidDynamicsApplication SHARED ${KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore)
target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore KratosSolidMechanicsApplication KratosMeshingApplication triangle tetgen)
##target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore KratosSolidMechanicsApplication triangle tetgen)
##target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore triangle tetgen)
##target_link_libraries(KratosPfemFluidDynamicsApplication KratosCore triangle tetgen msuite tetMeshOpt)
set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosPfemFluidDynamicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemFluidDynamicsApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemFluidDynamicsApplication.py" DESTINATION KratosMultiphysics )
