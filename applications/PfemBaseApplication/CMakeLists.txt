set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemBaseApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../PfemBaseApplication )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)


## generate variables with the sources
set( KRATOS_PFEM_BASE_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_base_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_base_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_base_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp

	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/modeler_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_data_transfer_utilities.cpp

	##conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/composite_condition.cpp

	##modelers
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/mesh_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/triangular_mesh_2D_modeler.cpp


)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemBaseApplication SHARED ${KRATOS_PFEM_BASE_APPLICATION_SOURCES})
#####***target_link_libraries(KratosPfemBaseApplication KratosCore)
##target_link_libraries(KratosPfemBaseApplication KratosCore KratosSolidMechanicsApplication KratosMeshingApplication triangle tetgen)
##target_link_libraries(KratosPfemBaseApplication KratosCore KratosSolidMechanicsApplication triangle tetgen)
target_link_libraries(KratosPfemBaseApplication KratosCore triangle tetgen)
##target_link_libraries(KratosPfemBaseApplication KratosCore triangle tetgen msuite tetMeshOpt)
set_target_properties(KratosPfemBaseApplication PROPERTIES PREFIX "")
install(TARGETS KratosPfemBaseApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemBaseApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemBaseApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemBaseApplication.py" DESTINATION KratosMultiphysics )
