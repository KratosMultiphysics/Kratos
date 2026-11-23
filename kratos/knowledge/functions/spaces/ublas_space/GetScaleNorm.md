---
type: C++ Function
title: GetScaleNorm
resource: spaces/ublas_space.h#L780-L800
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetDiagonalNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetMaxDiagonal
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double GetScaleNorm( const ProcessInfo& rProcessInfo, const MatrixType& rA, const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING )`

# Calls

- [GetDiagonalNorm](../../../functions/spaces/ublas_space/GetDiagonalNorm.md)
- [GetMaxDiagonal](../../../functions/spaces/ublas_space/GetMaxDiagonal.md)

# Called by

- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyConstraints.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyDirichletConditions.md)
- [CheckAndCorrectZeroDiagonalValues](../../../functions/spaces/ublas_space/CheckAndCorrectZeroDiagonalValues.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE.md)