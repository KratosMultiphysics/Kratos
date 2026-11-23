---
type: C++ Method
title: FixDof
resource: includes/dof.h#L334-L337
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/includes/node/Fix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/GenerateNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/ApplyDirichletConditions
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator
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

`void FixDof() noexcept`

# Called by

- [Fix](../../../../functions/includes/node/Fix.md)
- [GenerateNode](../../../../functions/modeler/structured_mesh_refinement_modeler/GenerateNode.md)
- [ApplyDirichletConditions](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/ApplyDirichletConditions.md)
- [MakePRestrictionOperator](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator.md)
- [BuildAndSolveLinearizedOnPreviousIteration](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/BuildAndSolveLinearizedOnPreviousIteration.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestScheme](../../../../functions/tests/cpp_tests/strategies/schemes/test_schemes/TestScheme.md)
- [ApplyFixity](../../../../functions/utilities/variable_utils/ApplyFixity.md)