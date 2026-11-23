---
type: C++ Function
title: CheckAndCorrectZeroDiagonalValues
resource: spaces/ublas_space.h#L715-L771
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spaces/ublas_space/GetScaleNorm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/parallel_utilities/IndexPartition/IndexPartition
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double CheckAndCorrectZeroDiagonalValues( const ProcessInfo& rProcessInfo, MatrixType& rA, VectorType& rb, const SCALING_DIAGONAL ScalingDiagonal = SCALING_DIAGONAL::NO_SCALING )`

# Calls

- [GetScaleNorm](../../../functions/spaces/ublas_space/GetScaleNorm.md)
- [IndexPartition](../../../functions/utilities/parallel_utilities/IndexPartition/IndexPartition.md)

# Called by

- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/ApplyDirichletConditions.md)
- [ApplyDirichletConditions](../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/ApplyDirichletConditions.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/spaces/test_ublas_space/KRATOS_TEST_CASE_IN_SUITE.md)