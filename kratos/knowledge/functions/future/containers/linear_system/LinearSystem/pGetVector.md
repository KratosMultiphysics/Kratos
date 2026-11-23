---
type: C++ Method
title: pGetVector
resource: future/containers/linear_system.h#L253-L259
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/containers/linear_system/LinearSystem/IsConsistent
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/amgcl_solver/AMGCLSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/iterative_solver/IterativeSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/PerformSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/python/add_containers_to_python/AddContainersToPython
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/ApplyLinearSystemConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateUpdateVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetResidualNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/EchoInfo
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/SolveSolutionStep
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/tests/cpp_tests/test_strategies/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`typename VectorType::Pointer pGetVector(DenseVectorTag Tag) const`

# Called by

- [IsConsistent](../../../../../functions/future/containers/linear_system/LinearSystem/IsConsistent.md)
- [PerformSolutionStep](../../../../../functions/future/linear_solvers/amgcl_solver/AMGCLSolver/PerformSolutionStep.md)
- [PerformSolutionStep](../../../../../functions/future/linear_solvers/iterative_solver/IterativeSolver/PerformSolutionStep.md)
- [PerformSolutionStep](../../../../../functions/future/linear_solvers/skyline_lu_factorization_solver/SkylineLUFactorizationSolver/PerformSolutionStep.md)
- [AddContainersToPython](../../../../../functions/future/python/add_containers_to_python/AddContainersToPython.md)
- [ApplyLinearSystemConstraints](../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyLinearSystemConstraints.md)
- [ApplyBlockBuildMasterSlaveConstraints](../../../../../functions/future/solving_strategies/builders/block_builder/BlockBuilder/ApplyBlockBuildMasterSlaveConstraints.md)
- [ApplyLinearSystemConstraints](../../../../../functions/future/solving_strategies/builders/elimination_builder/EliminationBuilder/ApplyLinearSystemConstraints.md)
- [CalculateUpdateVector](../../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/CalculateUpdateVector.md)
- [Update](../../../../../functions/future/solving_strategies/schemes/static_scheme/StaticScheme/Update.md)
- [GetResidualNorm](../../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/GetResidualNorm.md)
- [EchoInfo](../../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/EchoInfo.md)
- [SolveSolutionStep](../../../../../functions/future/solving_strategies/strategies/linear_strategy/LinearStrategy/SolveSolutionStep.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/future/tests/cpp_tests/test_schemes/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/future/tests/cpp_tests/test_strategies/KRATOS_TEST_CASE_IN_SUITE.md)