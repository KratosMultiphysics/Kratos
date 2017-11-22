set(CMAKE_INCLUDE_CURRENT_DIR ON)

## Kratos main source code
set( KRATOS_CORE_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/sources/global_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/deprecated_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/c2c_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/cfd_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/ale_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/mapping_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/dem_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/fsi_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/mat_variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/legacy_structural_app_vars.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/variables.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/integration_rules.cpp;
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
	${CMAKE_CURRENT_SOURCE_DIR}/processes/simple_mortar_mapper_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/containers/model.cpp;
	#${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_exception.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/containers/variable_data.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_local_smoothing_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_worst_element_smoothing_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_quality_weighted_smoothing_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/measure_mesh_quality_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/structured_mesh_generator_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_coarsening_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_node_collapsing_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/structured_mesh_generator_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_edge_swapping_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/find_intersected_geometrical_objects_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/calculate_discontinuous_distance_to_skin_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/calculate_distance_to_skin_process.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/processes/reorder_and_optimize_modelpart_process.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/utilities/exact_mortar_segmentation_utility.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sources/bounding_volume_tree.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/exception.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/code_location.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part_io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/reorder_consecutive_model_part_io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/io.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/sources/prime_numbers.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_geometry.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_triangle_2d_3.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_tetrahedra_3d_4.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modeler/tetrahedra_ball.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modeler/tetrahedra_edge_shell.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/modified_shape_functions.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/ausas_modified_shape_functions.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/triangle_2d_3_modified_shape_functions.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.cpp;
)

## Kratos python interface code
set( KRATOS_PYTHON_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_banded_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_part_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_compressed_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_constant_matrix_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_points_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_deprecated_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_c2c_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_cfd_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_ale_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_mapping_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_dem_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_fsi_variables_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_mat_variables_to_python.cpp
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
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_amgcl_solver_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_process_info_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_market_interface_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_serializer_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_table_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_search_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_kratos_parameters_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_testing_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/python/add_logger_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_to_python.cpp;
	${CMAKE_CURRENT_SOURCE_DIR}/python/kratos_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/sources/connectivity_preserve_modeler.cpp;
)

## Kratos testing soruces (testing engine, test are below)
file(GLOB_RECURSE KRATOS_CORE_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testing/*.cpp)

## Kratos I/o sources
file(GLOB_RECURSE KRATOS_CORE_INPUT_OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/input_output/*.cpp)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
	file(GLOB_RECURSE KRATOS_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Define library KratosCore to be included in all of the others
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES} ${KRATOS_CORE_TESTING_SOURCES} ${KRATOS_TESTING_SOURCES} ${KRATOS_CORE_INPUT_OUTPUT_SOURCES})
if(USE_COTIRE MATCHES ON)
    cotire(KratosCore)
endif(USE_COTIRE MATCHES ON)

target_link_libraries(KratosCore ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gidpost )
set_target_properties(KratosCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT,API")
install(TARGETS KratosCore DESTINATION libs )

## Define library Kratos which defines the basic python interface
add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(Kratos KratosCore)
if(USE_COTIRE MATCHES ON)
    cotire(Kratos)
endif(USE_COTIRE MATCHES ON)
set_target_properties(Kratos PROPERTIES PREFIX "")

## Changing the .dll suffix to .pyd (For Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(Kratos PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

## Changing the .dylib suffix to .so (For MacOS)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(Kratos PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS Kratos DESTINATION libs )
install(TARGETS KratosCore DESTINATION libs)
if(${INSTALL_PYTHON_FILES} MATCHES ON)
	get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
	install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION kratos/ FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
