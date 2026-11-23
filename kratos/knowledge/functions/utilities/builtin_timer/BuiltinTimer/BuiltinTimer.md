---
type: C++ Method
title: BuiltinTimer
resource: utilities/builtin_timer.h#L32-L34
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMassMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildDampingMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildLHS
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/BuildAndSolve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`BuiltinTimer() : mStartTime{std::chrono::steady_clock::now()}`

# Called by

- [Build](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/Build.md)
- [BuildWithSafeAssemble](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildWithSafeAssemble.md)
- [BuildWithThreadLocal](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithThreadLocal.md)
- [BuildWithLocalAllocation](../../../../functions/future/solving_strategies/schemes/implicit_scheme/BuildWithLocalAllocation.md)
- [BuildMassMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMassMatrix.md)
- [BuildDampingMatrix](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildDampingMatrix.md)
- [BuildMasterSlaveConstraints](../../../../functions/future/solving_strategies/schemes/implicit_scheme/ImplicitScheme/BuildMasterSlaveConstraints.md)
- [ApplyConstraintsForPeriodicConditions](../../../../functions/processes/apply_periodic_boundary_condition_process/ApplyPeriodicConditionProcess/ApplyConstraintsForPeriodicConditions.md)
- [Build](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/Build.md)
- [BuildLHS](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildLHS.md)
- [BuildAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolve.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [BuildRHSAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildRHSAndSolve.md)
- [Build](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/Build.md)
- [BuildAndSolve](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/BuildAndSolve.md)
- [Build](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/Build.md)
- [BuildAndSolveWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints.md)
- [BuildWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints.md)
- [BuildRHSNoDirichlet](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_parallel_utilities/KRATOS_TEST_CASE_IN_SUITE.md)