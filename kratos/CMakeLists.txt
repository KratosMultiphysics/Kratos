set(CMAKE_INCLUDE_CURRENT_DIR ON)

# Mpi
if(${MPI_NEEDED} MATCHES ON )
  add_subdirectory(mpi)
endif(${MPI_NEEDED} MATCHES ON )

## Kratos main source code
set( KRATOS_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/exact_mortar_segmentation_utility.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/global_variables.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/deprecated_variables.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/c2c_variables.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/cfd_variables.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/mesh_moving_variables.cpp;
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
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/process_info.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/simple_mortar_mapper_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/io.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/model_part.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/containers/model.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/standard_linear_solver_factory.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/standard_preconditioner_factory.cpp;
    #${CMAKE_CURRENT_SOURCE_DIR}/sources/kratos_exception.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/containers/variable_data.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/containers/data_value_container.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/fast_transfer_between_model_parts_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_local_smoothing_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_worst_element_smoothing_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_quality_weighted_smoothing_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/measure_mesh_quality_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/structured_mesh_generator_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_coarsening_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/mesh_node_collapsing_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/tetrahedra_mesh_edge_swapping_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/find_intersected_geometrical_objects_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/calculate_discontinuous_distance_to_skin_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/calculate_distance_to_skin_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/reorder_and_optimize_modelpart_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/compute_nodal_gradient_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/find_nodal_h_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/skin_detection_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/replace_elements_and_condition_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/bounding_volume_tree.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/exception.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/code_location.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/prime_numbers.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/connectivity_preserve_modeler.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/sources/parallel_environment.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_geometry.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_triangle_2d_3.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/divide_tetrahedra_3d_4.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/embedded_skin_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/geometrical_sensitivity_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/brute_force_point_locator.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/variable_utils.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/sub_model_parts_list_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/assign_unique_model_part_collection_tag_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/read_materials_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/variable_redistribution_utility.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/utilities/auxiliar_model_part_utilities.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/tetrahedra_ball.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modeler/tetrahedra_edge_shell.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/ausas_modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/triangle_2d_3_modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/triangle_2d_3_ausas_modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/tetrahedra_3d_4_modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/modified_shape_functions/tetrahedra_3d_4_ausas_modified_shape_functions.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/elements/mesh_element.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/conditions/mesh_condition.cpp;
)

## Kratos python interface code
set( KRATOS_PYTHON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/processes/assign_scalar_field_to_conditions_process.cpp;
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_part_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_points_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_deprecated_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_c2c_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_cfd_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_mesh_moving_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_mapping_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_dem_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_fsi_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_mat_variables_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_legacy_structural_app_vars_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_containers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_constitutive_law_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_geometries_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_kernel_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_vector_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_kratos_application_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_matrix_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_mesh_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_modeler_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_properties_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_node_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_amgcl_solver_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_linear_solvers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_factories_to_python.cpp
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
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_model_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_memory_info_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_response_functions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/kratos_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_constraint_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_data_communicator_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/python/add_parallel_environment_to_python.cpp
)

## Kratos testing engine soruces
file(GLOB_RECURSE KRATOS_CORE_TESTING_ENGINE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/testing/*.cpp)

## Kratos I/o sources
file(GLOB_RECURSE KRATOS_CORE_INPUT_OUTPUT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/input_output/*.cpp)

## Kratos tests sources. Enabled by default
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Define library KratosCore to be included in all of the others
add_library(KratosCore SHARED ${KRATOS_CORE_SOURCES} ${KRATOS_CORE_TESTING_ENGINE_SOURCES} ${KRATOS_CORE_INPUT_OUTPUT_SOURCES} ${KRATOS_TEST_SOURCES})
target_link_libraries(KratosCore PUBLIC gidpost) #${Boost_LIBRARIES} ${PYTHON_LIBRARIES} gidpost )
set_target_properties(KratosCore PROPERTIES COMPILE_DEFINITIONS "KRATOS_CORE=IMPORT,API")

## Define library Kratos which defines the basic python interface
pybind11_add_module(Kratos MODULE THIN_LTO ${KRATOS_PYTHON_SOURCES})
# add_library(Kratos SHARED ${KRATOS_PYTHON_SOURCES})
target_link_libraries(Kratos PRIVATE ${KRATOS_TEST_LIBRARIES} KratosCore)
# set_target_properties(Kratos PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    ## Kratos cannot be compiled with just one unity build as it takes way too much ram
    SET(COTIRE_MAXIMUM_NUMBER_OF_UNITY_INCLUDES "-j2")
    cotire(KratosCore)
    if(${KRATOS_BUILD_TESTING} MATCHES ON)
        cotire(KratosTests)
    endif(${KRATOS_BUILD_TESTING} MATCHES ON)

    ## Kratos cannot be compiled with just one unity build as it takes way too much ram
    SET(COTIRE_MAXIMUM_NUMBER_OF_UNITY_INCLUDES "-j2")

    ## This soruces in particular take too much ram by their own
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_amgcl_solver_to_python.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_strategies_to_python.cpp PROPERTIES COTIRE_EXCLUDED TRUE)
    set_source_files_properties (${CMAKE_CURRENT_SOURCE_DIR}/python/add_utilities_to_python.cpp PROPERTIES COTIRE_EXCLUDED TRUE)

    cotire(Kratos)
endif(USE_COTIRE MATCHES ON)

## Changing the .dll suffix to .pyd (For Windows)
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#     if(${KRATOS_BUILD_TESTING} MATCHES ON)
#         set_target_properties(KratosTests PROPERTIES SUFFIX .pyd)
#     endif(${KRATOS_BUILD_TESTING} MATCHES ON)
#     set_target_properties(Kratos PROPERTIES SUFFIX .pyd)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#
# ## Changing the .dylib suffix to .so (For MacOS)
# if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     if(${KRATOS_BUILD_TESTING} MATCHES ON)
#         set_target_properties(KratosTests PROPERTIES SUFFIX .so)
#     endif(${KRATOS_BUILD_TESTING} MATCHES ON)
#     set_target_properties(Kratos PROPERTIES SUFFIX .so)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS Kratos DESTINATION libs )
install(TARGETS KratosCore DESTINATION libs)

# Kratos Testing. Install everything except sources
if(${INSTALL_TESTING_FILES} MATCHES ON )
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION kratos
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
    )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install python scripts
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION kratos/ FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
