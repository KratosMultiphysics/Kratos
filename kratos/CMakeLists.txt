set(CMAKE_INCLUDE_CURRENT_DIR ON)


# ## add the include directories needed for the compilation
# message("current dir = " ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/sources )
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/container )
# include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/python )


## generate variables with the sources
set( KRATOS_CORE_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/sources/global_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/deprecated_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/c2c_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/cfd_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/dem_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/legacy_structural_app_vars.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/integration_rules.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/containers/variables_list_data_value_container.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/mmio.c;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/timer.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/serializer.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/kernel.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/periodic_condition.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/constitutive_law.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_components.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part_io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/reorder_consecutive_model_part_io.cpp;
	# ${CMAKE_CURRENT_SOURCE_DIR}/processes/calculate_signed_distance_to_3d_skin_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/utilities/logger.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_exception.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/containers/variable_data.cpp;
	#${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_smoothing_process.cpp;
	#${CMAKE_CURRENT_SOURCE_DIR}/processes/measure_mesh_quality_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/bounding_volume_tree.cpp;
     )

set( KRATOS_PYTHON_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_banded_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_part_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_compressed_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_constant_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_points_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_deprecated_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_c2c_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_cfd_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_dem_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_legacy_structural_app_vars_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_containers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_constitutive_law_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_coordinate_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_quadratures_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_geometries_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_scalar_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_hermitian_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_sparse_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_identity_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_integration_points_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_symmetric_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_io_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_triangular_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kernel_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_vector_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kratos_application_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_zero_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_mesh_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_modeler_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/vectorial_spaces_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_properties_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_node_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_process_info_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_parallel_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_parallel_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_market_interface_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_serializer_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_table_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_search_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kratos_parameters_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/kratos_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sources/mpi_connectivity_preserve_modeler.cpp;
	# ${CMAKE_CURRENT_SOURCE_DIR}/processes/mdpa_reorder_consecutive_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part_io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/reorder_consecutive_model_part_io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/io.cpp;
	#${CMAKE_CURRENT_SOURCE_DIR}/utilities/logger.cpp;
	# ${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_exception.cpp;

     )

###############################################################
## define library KratosCore to be included in all of the others
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES})
target_link_libraries(KratosCore ${Boost_PYTHON_LIBRARIES} ${PYTHON_LIBRARIES} gidpost )
set_target_properties(KratosCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT,API")
install(TARGETS KratosCore DESTINATION libs )

# if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#     SET_TARGET_PROPERTIES( KratosCore PROPERTIES COMPILE_FLAGS -fPIC)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION kratos/  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

###############################################################
## define library Kratos which defines the basic python interface
add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(Kratos KratosCore)
set_target_properties(Kratos PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(Kratos PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(Kratos PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS Kratos DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("subdir inc_dirs = ${inc_dirs}")
