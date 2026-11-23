---
type: C++ Method
title: SetBool
resource: sources/kratos_parameters.cpp#L842-L845
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetResetSolverIndexEachTry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/CheckAndAssignSettings
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/strategies/line_search_strategy/OverrideDefaultSettingsWithParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/__base_test_calculate_distance_to_path_process
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_2d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_3d
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofArrayUtilities
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/_SetUpTestModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_kratos_parameters/TestParameters/test_set_methods
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_mortar_mapper/TestMortarMapperCore/__base_test_mapping
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_obj_io/ReadModelPartFromOBJ
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Parameters::SetBool(const bool Value)`

# Called by

- [SetResetSolverIndexEachTry](../../../../functions/linear_solvers/fallback_linear_solver/FallbackLinearSolver/SetResetSolverIndexEachTry.md)
- [ImportModelPart](../../../../functions/mpi/python_scripts/distributed_import_model_part_utility/DistributedImportModelPartUtility/ImportModelPart.md)
- [__execute_test](../../../../functions/mpi/tests/test_distributed_import_model_part_utility/TestDistributedImportModelPartUtility/__execute_test.md)
- [CheckAndAssignSettings](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/CheckAndAssignSettings.md)
- [SimpleMortarMapperProcessWrapper](../../../../functions/processes/simple_mortar_mapper_wrapper_process/SimpleMortarMapperProcessWrapper/SimpleMortarMapperProcessWrapper.md)
- [__init__](../../../../functions/python_scripts/assign_vector_variable_process/AssignVectorVariableProcess/__init__.md)
- [ReplaceDeprecatedVariableName](../../../../functions/python_scripts/deprecation_management/DeprecationManager/ReplaceDeprecatedVariableName.md)
- [__init__](../../../../functions/python_scripts/vtk_output_process/VtkOutputProcess/__init__.md)
- [OverrideDefaultSettingsWithParameters](../../../../functions/solving_strategies/strategies/line_search_strategy/OverrideDefaultSettingsWithParameters.md)
- [AddBool](../../../../functions/sources/kratos_parameters/Parameters/AddBool.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_kratos_parameters/KRATOS_TEST_CASE_IN_SUITE.md)
- [__base_test_calculate_distance_to_path_process](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/__base_test_calculate_distance_to_path_process.md)
- [test_naca_0012_calculate_distance_to_skin_2d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_2d.md)
- [test_naca_0012_calculate_distance_to_skin_3d](../../../../functions/tests/test_calculate_distance_to_skin/TestCalculateDistanceToSkin/test_naca_0012_calculate_distance_to_skin_3d.md)
- [testDofArrayUtilities](../../../../functions/tests/test_dofs/TestDofs/testDofArrayUtilities.md)
- [_SetUpTestModelPart](../../../../functions/tests/test_edge_based_gradient_recovery_process/TestEdgeBasedGradientRecoveryProcess/_SetUpTestModelPart.md)
- [test_embedded_skin_mapping](../../../../functions/tests/test_embedded_skin_mapping/TestEmbeddedSkinMapping/test_embedded_skin_mapping.md)
- [test_set_methods](../../../../functions/tests/test_kratos_parameters/TestParameters/test_set_methods.md)
- [setUp](../../../../functions/tests/test_mls_shape_functions_utility/TestMLSShapeFunctionsUtility/setUp.md)
- [__base_test_mapping](../../../../functions/tests/test_mortar_mapper/TestMortarMapperCore/__base_test_mapping.md)
- [ReadModelPartFromOBJ](../../../../functions/tests/test_obj_io/ReadModelPartFromOBJ.md)
- [setUp](../../../../functions/tests/test_rbf_shape_functions_utility/TestRBFShapeFunctionsUtility/setUp.md)