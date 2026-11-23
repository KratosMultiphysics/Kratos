---
type: C++ Method
title: FreeDof
resource: includes/dof.h#L341-L344
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/node/Free
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/schemes/test_schemes/TestScheme
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/ApplyFixity
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FreeDof() noexcept`

# Called by

- [Free](../../../../functions/includes/node/Free.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestScheme](../../../../functions/tests/cpp_tests/strategies/schemes/test_schemes/TestScheme.md)
- [ApplyFixity](../../../../functions/utilities/variable_utils/ApplyFixity.md)