---
type: C++ Method
title: SetEquationId
resource: includes/dof.h#L320-L323
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/mpi/includes/mpi_communicator/SendTools< Node/DofsContainerType >/ReadBuffer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSetEquationIds
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Finalize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetUpSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/SetUpSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpSystemWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/linear_solvers/test_monotonicity_solver/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/convergence_criteria/test_displacement_criteria/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetEquationId(EquationIdType NewEquationId) noexcept`

# Called by

- [ReadBuffer](../../../../functions/mpi/includes/mpi_communicator/SendTools< Node/DofsContainerType >/ReadBuffer.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/mpi/tests/cpp_tests/sources/test_mpi_communicator/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetUpDofSetEquationIds](../../../../functions/solving_strategies/builder_and_solvers/explicit_builder/ExplicitBuilder/SetUpDofSetEquationIds.md)
- [Initialize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Initialize.md)
- [Finalize](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/Finalize.md)
- [SetUpSystem](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/SetUpSystem.md)
- [MakePRestrictionOperator](../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_utilities/MakePRestrictionOperator.md)
- [SetUpSystem](../../../../functions/solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver/ResidualBasedBlockBuilderAndSolver/SetUpSystem.md)
- [SetUpSystem](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/SetUpSystem.md)
- [SetUpSystemWithConstraints](../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/SetUpSystemWithConstraints.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/linear_solvers/test_monotonicity_solver/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/sources/test_dof/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/strategies/convergence_criteria/test_displacement_criteria/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/strategies/convergence_criteria/test_residual_criteria/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_constraint_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [SetDofEquationIds](../../../../functions/utilities/dof_utilities/dof_array_utilities/DofArrayUtilities/SetDofEquationIds.md)