---
type: C++ Method
title: CleanFunctionName
resource: sources/code_location.cpp#L59-L83
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/code_location/CodeLocation/RemoveNamespace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/ReduceTemplateArgumentsToFirstN
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/code_location/CodeLocation/ReplaceAll
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeIndependentResidual
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint/Create
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildRHSAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/ApplyRHSConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/BuildLHS_CompleteOnFreeRows
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/BuildAndSolveLinearizedOnPreviousIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_code_location/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/profiler/Profiler/Write
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string CodeLocation::CleanFunctionName() const`

# Calls

- [RemoveNamespace](../../../../functions/sources/code_location/CodeLocation/RemoveNamespace.md)
- [ReduceTemplateArgumentsToFirstN](../../../../functions/sources/code_location/CodeLocation/ReduceTemplateArgumentsToFirstN.md)
- [ReplaceAll](../../../../functions/sources/code_location/CodeLocation/ReplaceAll.md)

# Called by

- [Create](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/linear_multifreedom_constraint/Create.md)
- [ComputeIndependentResidual](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/ComputeIndependentResidual.md)
- [Create](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/multifreedom_constraint/Create.md)
- [BuildRHSAndSolve](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/BuildRHSAndSolve.md)
- [ApplyRHSConstraints](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/ApplyRHSConstraints.md)
- [BuildLHS_CompleteOnFreeRows](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/BuildLHS_CompleteOnFreeRows.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_code_location/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_exceptions/KRATOS_TEST_CASE_IN_SUITE.md)
- [Write](../../../../functions/utilities/profiler/Profiler/Write.md)