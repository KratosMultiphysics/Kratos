---
type: C++ Method
title: ValidateAndAssignDefaults
resource: sources/kratos_parameters.cpp#L1345-L1353
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/AddMissingParameters
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/linear_solvers/amgcl_solver/AMGCLSolver/AMGCLSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/cg_solver/CGSolver/CGSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/linear_solver/LinearSolver/LinearSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/SkylineLUFactorizationSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/BlockBuilder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/builder/Builder/Builder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/EliminationBuilder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/ImplicitStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/LinearStrategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/strategy/Strategy/Strategy
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/logger_table_output/LoggerTableOutput/LoggerTableOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/obj_io/ObjIO/ObjIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/stl_io/StlIO/StlIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtk_output/VtkOutput/VtkOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/AMGCL_NS_Solver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/DeflatedCGSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/iterative_solver/IterativeSolver/IterativeSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/monotonicity_preserving_solver/MonotonicityPreservingSolver/MonotonicityPreservingSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/power_iteration_eigenvalue_solver/PowerIterationEigenvalueSolver/PowerIterationEigenvalueSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/power_iteration_highest_eigenvalue_solver/PowerIterationHighestEigenvalueSolver/PowerIterationHighestEigenvalueSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/RayleighQuotientIterationEigenvalueSolver
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/coloring/voxel_mesher_coloring/VoxelMesherColoring/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CombineModelPartModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ConnectivityPreserveModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation/VoxelMesherKeyPlaneGeneration/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_by_bounding_box/KeyPlaneGenerationByBoundingBox/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_by_outer_shell/KeyPlaneGenerationByOuterShell/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/operation/voxel_mesher_operation/VoxelMesherOperation/ValidateParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/VoxelMeshGeneratorModeler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/ApplyConstantScalarValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/ApplyRayCastingInterfaceRecognitionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ApplyRayCastingProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/AssignScalarFieldToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess<TEntity, THistorical>/AssignScalarVariableToEntitiesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceToPathProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToSkinProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/CalculateNodalDistanceToSkinProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/FindConditionsNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/FindGlobalNodalEntityNeighboursProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/FluxCorrectedTransportConvectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FromJSONCheckResultProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/IntegrationValuesExtrapolationToNodesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ParallelDistanceCalculationProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SkinDetectionProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/TetrahedralMeshOrientationCheck
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_flag_process/AssignFlagProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_conditions_process/AssignScalarInputToConditionsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_elements_process/AssignScalarInputToElementsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_input_to_nodes_process/AssignScalarInputToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_conditions_process/AssignScalarVariableToConditionsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_constraints_process/AssignScalarVariableToConstraintsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_elements_process/AssignScalarVariableToElementsProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_scalar_variable_to_nodes_process/AssignScalarVariableToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_condition_process/AssignVectorByDirectionToConditionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_element_process/AssignVectorByDirectionToElementProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_by_direction_to_node_process/AssignVectorByDirectionToNodeProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_conditions_process/AssignVectorVariableToConditionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_constraints_process/AssignVectorVariableToConstraintProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_elements_process/AssignVectorVariableToElementProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_to_nodes_process/AssignVectorVariableToNodeProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_scalar_to_nodes_process/CheckScalarToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/csv_points_output_process/Factory
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/json_output_process/JsonOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/line_output_process/LineOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/output_process_with_controller/OutputProcessWithController/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/python_solver/PythonSolver/ValidateSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/read_csv_table_utility/ReadCsvTableUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/restart_utility/RestartUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/save_restart_process/SaveRestartProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/skin_detection_process/SkinDetectionProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/timer_process/TimerProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/builder_and_solver/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/AssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResidualBasedEliminationBuilderAndSolverComponentwise
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/ResidualBasedAdjointBossakScheme
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incremental_aitken_static_scheme/ResidualBasedIncrementalAitkenStaticScheme/ResidualBasedIncrementalAitkenStaticScheme
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/ResidualBasedIncrementalUpdateStaticSchemeSlip
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/scheme/Scheme/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/ValidateAndAssignParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_validation_fails_due_to_wrong_type
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_validation_fails_due_to_wrong_spelling
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_validation_succeds_error_on_first_level
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_validation_succeeds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_null_vs_null_validation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_double_vs_null_validation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/ShiftedBoundaryMeshlessInterfaceUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::ValidateAndAssignDefaults(const Parameters& rDefaultParameters)`

# Calls

- [AddMissingParameters](../../../../functions/sources/kratos_parameters/Parameters/AddMissingParameters.md)

# Called by

- [AMGCLSolver](../../../../functions/future/linear_solvers/amgcl_solver/AMGCLSolver/AMGCLSolver.md)
- [CGSolver](../../../../functions/future/linear_solvers/cg_solver/CGSolver/CGSolver.md)
- [LinearSolver](../../../../functions/future/linear_solvers/linear_solver/LinearSolver/LinearSolver.md)
- [SkylineLUFactorizationSolver](../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/SkylineLUFactorizationSolver.md)
- [BlockBuilder](../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/BlockBuilder.md)
- [Builder](../../../../functions/future/solving_strategies/builders/builder/Builder/Builder.md)
- [EliminationBuilder](../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/EliminationBuilder.md)
- [ValidateAndAssignParameters](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/ValidateAndAssignParameters.md)
- [ImplicitStrategy](../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/ImplicitStrategy.md)
- [LinearStrategy](../../../../functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/LinearStrategy.md)
- [Strategy](../../../../functions/future/solving_strategies/strategies/strategy/Strategy/Strategy.md)
- [LoggerTableOutput](../../../../functions/input_output/logger_table_output/LoggerTableOutput/LoggerTableOutput.md)
- [ObjIO](../../../../functions/input_output/obj_io/ObjIO/ObjIO.md)
- [StlIO](../../../../functions/input_output/stl_io/StlIO/StlIO.md)
- [VtkOutput](../../../../functions/input_output/vtk_output/VtkOutput/VtkOutput.md)
- [AMGCL_NS_Solver](../../../../functions/linear_solvers/amgcl_ns_solver/AMGCL_NS_Solver/AMGCL_NS_Solver.md)
- [ApplySettings](../../../../functions/linear_solvers/amgcl_solver_impl/AMGCLSolver<TSparse,TDense>/ApplySettings.md)
- [DeflatedCGSolver](../../../../functions/linear_solvers/deflated_cg_solver/DeflatedCGSolver/DeflatedCGSolver.md)
- [FallbackLinearSolver](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/FallbackLinearSolver.md)
- [IterativeSolver](../../../../functions/linear_solvers/iterative_solver/IterativeSolver/IterativeSolver.md)
- [MonotonicityPreservingSolver](../../../../functions/linear_solvers/monotonicity_preserving_solver/MonotonicityPreservingSolver/MonotonicityPreservingSolver.md)
- [PowerIterationEigenvalueSolver](../../../../functions/linear_solvers/power_iteration_eigenvalue_solver/PowerIterationEigenvalueSolver/PowerIterationEigenvalueSolver.md)
- [PowerIterationHighestEigenvalueSolver](../../../../functions/linear_solvers/power_iteration_highest_eigenvalue_solver/PowerIterationHighestEigenvalueSolver/PowerIterationHighestEigenvalueSolver.md)
- [RayleighQuotientIterationEigenvalueSolver](../../../../functions/linear_solvers/rayleigh_quotient_iteration_eigenvalue_solver/RayleighQuotientIterationEigenvalueSolver/RayleighQuotientIterationEigenvalueSolver.md)
- [CleanUpProblematicTrianglesModeler](../../../../functions/modeler/clean_up_problematic_triangles_modeler/CleanUpProblematicTrianglesModeler.md)
- [ValidateParameters](../../../../functions/modeler/coloring/voxel_mesher_coloring/VoxelMesherColoring/ValidateParameters.md)
- [CombineModelPartModeler](../../../../functions/modeler/combine_model_part_modeler/CombineModelPartModeler/CombineModelPartModeler.md)
- [ConnectivityPreserveModeler](../../../../functions/modeler/connectivity_preserve_modeler/ConnectivityPreserveModeler/ConnectivityPreserveModeler.md)
- [CopyPropertiesModeler](../../../../functions/modeler/copy_properties_modeler/CopyPropertiesModeler/CopyPropertiesModeler.md)
- [CreateEntitiesFromGeometriesModeler](../../../../functions/modeler/create_entities_from_geometries_modeler/CreateEntitiesFromGeometriesModeler.md)
- [ValidateParameters](../../../../functions/modeler/entity_generation/voxel_mesher_entity_generation/VoxelMesherEntityGeneration/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/key_plane_generation/key_plane_generation/VoxelMesherKeyPlaneGeneration/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/key_plane_generation/key_plane_generation_by_bounding_box/KeyPlaneGenerationByBoundingBox/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/key_plane_generation/key_plane_generation_by_outer_shell/KeyPlaneGenerationByOuterShell/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/key_plane_generation/key_plane_generation_with_refinement/KeyPlaneGenerationWithRefinement/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/operation/compute_surrogate_boundary_data/ComputeSurrogateBoundaryData/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/operation/find_contacts_in_skin_model_part/FindContactsInSkinModelPart/ValidateParameters.md)
- [ValidateParameters](../../../../functions/modeler/operation/voxel_mesher_operation/VoxelMesherOperation/ValidateParameters.md)
- [VoxelMeshGeneratorModeler](../../../../functions/modeler/voxel_mesh_generator_modeler/VoxelMeshGeneratorModeler/VoxelMeshGeneratorModeler.md)
- [ApplyConstantScalarValueProcess](../../../../functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/ApplyConstantScalarValueProcess.md)
- [ApplyConstantVectorValueProcess](../../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/ApplyConstantVectorValueProcess.md)
- [ApplyRayCastingInterfaceRecognitionProcess](../../../../functions/processes/apply_ray_casting_interface_recognition_process/ApplyRayCastingInterfaceRecognitionProcess<TDim>/ApplyRayCastingInterfaceRecognitionProcess.md)
- [ApplyRayCastingProcess](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/ApplyRayCastingProcess.md)
- [AssignScalarFieldToEntitiesProcess](../../../../functions/processes/assign_scalar_field_to_entities_process/AssignScalarFieldToEntitiesProcess<TEntity, THistorical>/AssignScalarFieldToEntitiesProcess.md)
- [AssignScalarInputToEntitiesProcess](../../../../functions/processes/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess<TEntity, THistorical>/AssignScalarInputToEntitiesProcess.md)
- [AssignScalarVariableToEntitiesProcess](../../../../functions/processes/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess<TEntity, THistorical>/AssignScalarVariableToEntitiesProcess.md)
- [CalculateDistanceToPathProcess](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceToPathProcess.md)
- [CalculateDistanceToSkinProcess](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToSkinProcess.md)
- [CalculateEmbeddedNodalVariableFromSkinProcess](../../../../functions/processes/calculate_embedded_nodal_variable_from_skin_process/CalculateEmbeddedNodalVariableFromSkinProcess/CalculateEmbeddedNodalVariableFromSkinProcess.md)
- [CalculateNodalDistanceToSkinProcess](../../../../functions/processes/calculate_nodal_distance_to_skin_process/CalculateNodalDistanceToSkinProcess/CalculateNodalDistanceToSkinProcess.md)
- [EdgeBasedGradientRecoveryProcess](../../../../functions/processes/edge_based_gradient_recovery_process/EdgeBasedGradientRecoveryProcess/EdgeBasedGradientRecoveryProcess.md)
- [FindConditionsNeighboursProcess](../../../../functions/processes/find_conditions_neighbours_process/FindConditionsNeighboursProcess/FindConditionsNeighboursProcess.md)
- [FindGlobalNodalEntityNeighboursProcess](../../../../functions/processes/find_global_nodal_entity_neighbours_process/FindGlobalNodalEntityNeighboursProcess<TContainerType>/FindGlobalNodalEntityNeighboursProcess.md)
- [FluxCorrectedTransportConvectionProcess](../../../../functions/processes/flux_corrected_transport_convection_process/FluxCorrectedTransportConvectionProcess/FluxCorrectedTransportConvectionProcess.md)
- [FromJSONCheckResultProcess](../../../../functions/processes/from_json_check_result_process/FromJSONCheckResultProcess/FromJSONCheckResultProcess.md)
- [IntegrationValuesExtrapolationToNodesProcess](../../../../functions/processes/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/IntegrationValuesExtrapolationToNodesProcess.md)
- [LevelSetConvectionProcess](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/LevelSetConvectionProcess.md)
- [ParallelDistanceCalculationProcess](../../../../functions/processes/parallel_distance_calculation_process/ParallelDistanceCalculationProcess<TDim>/ParallelDistanceCalculationProcess.md)
- [SimpleMortarMapperProcessWrapper](../../../../functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper.md)
- [SkinDetectionProcess](../../../../functions/processes/skin_detection_process/SkinDetectionProcess<TDim>/SkinDetectionProcess.md)
- [SplitInternalInterfacesProcess](../../../../functions/processes/split_internal_interfaces_process/SplitInternalInterfacesProcess.md)
- [ValidateTheDefaultParameters](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/ValidateTheDefaultParameters.md)
- [TetrahedralMeshOrientationCheck](../../../../functions/processes/tetrahedral_mesh_orientation_check/TetrahedralMeshOrientationCheck/TetrahedralMeshOrientationCheck.md)
- [VariationalDistanceCalculationProcess](../../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/VariationalDistanceCalculationProcess.md)
- [__init__](../../../../functions/python_scripts/apply_periodic_condition_process/ApplyPeriodicConditionProcessPythonWrapper/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_flag_process/AssignFlagProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_master_slave_constraints_to_neighbours_process/AssignMasterSlaveConstraintsToNeighboursProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_conditions_process/AssignScalarInputToConditionsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_elements_process/AssignScalarInputToElementsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_entities_process/AssignScalarInputToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_input_to_nodes_process/AssignScalarInputToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_process/AssignScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_conditions_process/AssignScalarVariableToConditionsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_constraints_process/AssignScalarVariableToConstraintsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_elements_process/AssignScalarVariableToElementsProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_entities_process/AssignScalarVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_scalar_variable_to_nodes_process/AssignScalarVariableToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_process/AssignVectorByDirectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_condition_process/AssignVectorByDirectionToConditionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_element_process/AssignVectorByDirectionToElementProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_entity_process/AssignVectorByDirectionToEntityProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_by_direction_to_node_process/AssignVectorByDirectionToNodeProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_conditions_process/AssignVectorVariableToConditionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_constraints_process/AssignVectorVariableToConstraintProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_elements_process/AssignVectorVariableToElementProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_entities_process/AssignVectorVariableToEntitiesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_to_nodes_process/AssignVectorVariableToNodeProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_scalar_base_process/CheckScalarBaseProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_scalar_to_nodes_process/CheckScalarToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/check_vector_components_to_nodes_process/CheckVectorComponentsToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/compare_two_files_check_process/CompareTwoFilesCheckProcess/__init__.md)
- [Factory](../../../../functions/python_scripts/csv_points_output_process/Factory.md)
- [__init__](../../../../functions/python_scripts/fix_scalar_variable_process/FixScalarVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/fix_vector_variable_process/FixVectorVariableProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_hdf5_postprocess/GiDHDF5OutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/__init__.md)
- [ExecuteInitialize](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)
- [_initialize_cut_output](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)
- [__init__](../../../../functions/python_scripts/integration_values_extrapolation_to_nodes_process/IntegrationValuesExtrapolationToNodesProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/json_output_process/JsonOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/line_output_process/LineOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/multiple_points_output_process/MultiplePointsOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/output_process_with_controller/OutputProcessWithController/__init__.md)
- [__init__](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/print_info_in_file_process/PrintInfoInFileProcess/__init__.md)
- [ValidateSettings](../../../../functions/python_scripts/python_solver/PythonSolver/ValidateSettings.md)
- [__init__](../../../../functions/python_scripts/read_csv_table_utility/ReadCsvTableUtility/__init__.md)
- [__init__](../../../../functions/python_scripts/restart_utility/RestartUtility/__init__.md)
- [__init__](../../../../functions/python_scripts/save_restart_process/SaveRestartProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/skin_detection_process/SkinDetectionProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/sub_model_part_entities_boolean_operation_process/SubModelPartEntitiesBooleanOperationProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/tikz_output_process/TikZOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/time_based_ascii_file_writer_utility/TimeBasedAsciiFileWriterUtility/__init__.md)
- [__init__](../../../../functions/python_scripts/timer_process/TimerProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/unv_output_process/UnvOutputProcess/__init__.md)
- [__init__](../../../../functions/python_scripts/vtu_output_process/VtuOutputProcess/__init__.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/ValidateAndAssignParameters.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/ValidateAndAssignParameters.md)
- [AugmentedLagrangeConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/AugmentedLagrangeConstraintAssembler.md)
- [MasterSlaveConstraintAssembler](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/MasterSlaveConstraintAssembler.md)
- [AssignSettings](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/AssignSettings.md)
- [ResidualBasedEliminationBuilderAndSolverComponentwise](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ResidualBasedEliminationBuilderAndSolverComponentwise.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/convergencecriterias/convergence_criteria/ConvergenceCriteria/ValidateAndAssignParameters.md)
- [ResidualBasedAdjointBossakScheme](../../../../functions/solving_strategies/schemes/residual_based_adjoint_bossak_scheme/ResidualBasedAdjointBossakScheme/ResidualBasedAdjointBossakScheme.md)
- [ResidualBasedIncrementalAitkenStaticScheme](../../../../functions/solving_strategies/schemes/residualbased_incremental_aitken_static_scheme/ResidualBasedIncrementalAitkenStaticScheme/ResidualBasedIncrementalAitkenStaticScheme.md)
- [ResidualBasedIncrementalUpdateStaticSchemeSlip](../../../../functions/solving_strategies/schemes/residualbased_incrementalupdate_static_scheme_slip/ResidualBasedIncrementalUpdateStaticSchemeSlip/ResidualBasedIncrementalUpdateStaticSchemeSlip.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/schemes/scheme/Scheme/ValidateAndAssignParameters.md)
- [ValidateAndAssignParameters](../../../../functions/solving_strategies/strategies/solving_strategy/SolvingStrategy/ValidateAndAssignParameters.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_validation_fails_due_to_wrong_type](../../../../functions/tests/test_kratos_parameters/TestParameters/test_validation_fails_due_to_wrong_type.md)
- [test_validation_fails_due_to_wrong_spelling](../../../../functions/tests/test_kratos_parameters/TestParameters/test_validation_fails_due_to_wrong_spelling.md)
- [test_validation_succeds_error_on_first_level](../../../../functions/tests/test_kratos_parameters/TestParameters/test_validation_succeds_error_on_first_level.md)
- [test_validation_succeeds](../../../../functions/tests/test_kratos_parameters/TestParameters/test_validation_succeeds.md)
- [test_null_vs_null_validation](../../../../functions/tests/test_kratos_parameters/TestParameters/test_null_vs_null_validation.md)
- [test_double_vs_null_validation](../../../../functions/tests/test_kratos_parameters/TestParameters/test_double_vs_null_validation.md)
- [SensitivityBuilder](../../../../functions/utilities/sensitivity_builder/SensitivityBuilder/SensitivityBuilder.md)
- [ShiftedBoundaryMeshlessInterfaceUtility](../../../../functions/utilities/shifted_boundary_meshless_interface_utility/ShiftedBoundaryMeshlessInterfaceUtility/ShiftedBoundaryMeshlessInterfaceUtility.md)