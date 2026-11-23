---
type: C++ Function
title: GetAveragevalueDiagonal
resource: spaces/ublas_space.h#L832-L835
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetMaxDiagonal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/GetMinDiagonal
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double GetAveragevalueDiagonal(const MatrixType& rA)`

# Calls

- [GetMaxDiagonal](../../../functions/spaces/ublas_space/GetMaxDiagonal.md)
- [GetMinDiagonal](../../../functions/spaces/ublas_space/GetMinDiagonal.md)

# Called by

- [ApplyConstraints](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver_with_lagrange_multiplier/ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier/ApplyConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE.md)