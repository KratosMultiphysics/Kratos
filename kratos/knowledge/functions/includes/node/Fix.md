---
type: C++ Function
title: Fix
resource: includes/node.h#L851-L866
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/dof/Dof/FixDof
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/openmp_utils/IsInParallel
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/tests/cpp_tests/test_strategies/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/json_io/CheckAndAssignNodalValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/json_io/CheckAndAssignNodalComponentValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/utilities/test_debug_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/InternalApplyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/InternalApplyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/add_node_to_python/NodeFix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/read_modelpart_from_json/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreFixityData
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_flux_corrected_transport_convection_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementAllDoFsMaster
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementWithZeroContribution
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/ExtendedTestBuilderAndSolverDisplacement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/GenerateTestExplicitBuilderModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_dofs/TestDofs/testDofExport
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_BFECC
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_BFECC_algebraic
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_tau_nodal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_linear_constraints/TestLinearConstraints/test_constraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_tensor_adaptors/fixity_setter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ApplyMinimalScalarFixity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void Fix(const TVariableType& rDofVariable)`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [FixDof](../../../functions/includes/dof/Dof/FixDof.md)
- [IsInParallel](../../../functions/utilities/openmp_utils/IsInParallel.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/future/tests/cpp_tests/test_strategies/KRATOS_TEST_CASE_IN_SUITE.md)
- [CheckAndAssignNodalValue](../../../functions/includes/json_io/CheckAndAssignNodalValue.md)
- [CheckAndAssignNodalComponentValue](../../../functions/includes/json_io/CheckAndAssignNodalComponentValue.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/mpi/tests/cpp_tests/utilities/test_debug_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [InternalApplyValue](../../../functions/processes/apply_constant_scalarvalue_process/ApplyConstantScalarValueProcess/InternalApplyValue.md)
- [InternalApplyValue](../../../functions/processes/apply_constant_vectorvalue_process/ApplyConstantVectorValueProcess/InternalApplyValue.md)
- [Execute](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/Execute.md)
- [SynchronizeFixity](../../../functions/processes/variational_distance_calculation_process/VariationalDistanceCalculationProcess/SynchronizeFixity.md)
- [NodeFix](../../../functions/python/add_node_to_python/NodeFix.md)
- [ReadModelPart](../../../functions/python_scripts/read_modelpart_from_json/ReadModelPart.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/InitializeSolutionStep.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bdf_displacement_scheme/ResidualBasedBDFDisplacementScheme/InitializeSolutionStep.md)
- [InitializeSolutionStep](../../../functions/solving_strategies/schemes/residual_based_bossak_displacement_scheme/ResidualBasedBossakDisplacementScheme/InitializeSolutionStep.md)
- [ReadNodalDofVariableData](../../../functions/sources/model_part_io/ModelPartIO/ReadNodalDofVariableData.md)
- [StoreFixityData](../../../functions/tensor_adaptors/fixity_tensor_adaptor/FixityTensorAdaptor/StoreFixityData.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/processes/test_flux_corrected_transport_convection_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE.md)
- [BasicTestBuilderAndSolverDisplacement](../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacement.md)
- [BasicTestBuilderAndSolverDisplacementAllDoFsMaster](../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementAllDoFsMaster.md)
- [BasicTestBuilderAndSolverDisplacementWithZeroContribution](../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/BasicTestBuilderAndSolverDisplacementWithZeroContribution.md)
- [ExtendedTestBuilderAndSolverDisplacement](../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_builder_and_solver/ExtendedTestBuilderAndSolverDisplacement.md)
- [GenerateTestExplicitBuilderModelPart](../../../functions/tests/cpp_tests/strategies/builder_and_solvers/test_explicit_builder_and_solver/GenerateTestExplicitBuilderModelPart.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [testDofExport](../../../functions/tests/test_dofs/TestDofs/testDofExport.md)
- [test_graph_and_connected_components_non_consecutive](../../../functions/tests/test_graph_utilities/TestGraphUtilities/test_graph_and_connected_components_non_consecutive.md)
- [test_levelset_convection](../../../functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection.md)
- [test_levelset_convection_BFECC](../../../functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_BFECC.md)
- [test_levelset_convection_BFECC_algebraic](../../../functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_BFECC_algebraic.md)
- [test_levelset_convection_tau_nodal](../../../functions/tests/test_levelset_convection/TestLevelSetConvection/test_levelset_convection_tau_nodal.md)
- [test_constraints](../../../functions/tests/test_linear_constraints/TestLinearConstraints/test_constraints.md)
- [fixity_setter](../../../functions/tests/test_tensor_adaptors/fixity_setter.md)
- [ApplyMinimalScalarFixity](../../../functions/utilities/model_part_graph_utilities/ModelPartGraphUtilities/ApplyMinimalScalarFixity.md)