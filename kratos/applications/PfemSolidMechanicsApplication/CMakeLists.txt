set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemSolidMechanicsApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../SolidMechanicsApplication )
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
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_solid_mechanics_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/skin_multiple_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/wall_tip_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/contact_domain_2D_condition.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/axisym_contact_domain_2D_condition.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/modeler_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_data_transfer_utilities.cpp

	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/triangle_mesh_2D_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/contact_domain_2D_modeler.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemSolidMechanicsApplication SHARED ${KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore)
target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore KratosSolidMechanicsApplication KratosMeshingApplication triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore KratosSolidMechanicsApplication triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore triangle tetgen)
##target_link_libraries(KratosPfemSolidMechanicsApplication KratosCore triangle tetgen msuite tetMeshOpt)
set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES PREFIX "")
install(TARGETS KratosPfemSolidMechanicsApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemSolidMechanicsApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemSolidMechanicsApplication.py" DESTINATION KratosMultiphysics )
