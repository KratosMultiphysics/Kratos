---
type: C++ Method
title: HasReaction
resource: includes/dof.h#L363-L366
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/variables_list/pGetDofReaction
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/Check
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/SetUpDofSet
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSetWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool HasReaction() const`

# Calls

- [pGetDofReaction](../../../../functions/containers/variables_list/pGetDofReaction.md)

# Called by

- [Check](../../../../functions/future/solving_strategies/strategies/implicit_strategy/ImplicitStrategy/Check.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSet.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpDofSet.md)
- [SetUpDofSet](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/SetUpDofSet.md)
- [SetUpDofSetWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpDofSetWithConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_node/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetUpDofArray](../../../../functions/utilities/dof_utilities/block_build_dof_array_utility/BlockBuildDofArrayUtility/SetUpDofArray.md)
- [SetUpDofArray](../../../../functions/utilities/dof_utilities/elimination_build_dof_array_utility/EliminationBuildDofArrayUtility/SetUpDofArray.md)