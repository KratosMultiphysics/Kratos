---
type: C++ Method
title: Tag
resource: solving_strategies/builder_and_solvers/p_multigrid/status_stream.hpp#L43-L46
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler/NoOpConstraintAssembler/FinalizeConstraintIteration
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ExecuteConstraintLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::pair<Report,int> Tag(int ReportVerbosity)`

# Called by

- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler/NoOpConstraintAssembler/FinalizeConstraintIteration.md)
- [ExecuteConstraintLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ExecuteConstraintLoop.md)
- [ExecuteMultigridLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop.md)
- [ExecuteConstraintLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop.md)