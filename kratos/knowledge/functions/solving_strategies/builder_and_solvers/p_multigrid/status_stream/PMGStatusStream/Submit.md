---
type: C++ Method
title: Submit
resource: solving_strategies/builder_and_solvers/p_multigrid/status_stream.cpp#L314-L323
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/DuplicateFilter
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState
    resolved_by: tree-sitter
    confidence: exact
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

`void PMGStatusStream::Submit(const TaggedReport& rReport, int SubmitterVerbosity)`

# Calls

- [DuplicateFilter](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/DuplicateFilter.md)
- [WriteIntermediateState](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/status_stream/PMGStatusStream/Impl/WriteIntermediateState.md)

# Called by

- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/augmented_lagrange_constraint_assembler/AugmentedLagrangeConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/master_slave_constraint_assembler/MasterSlaveConstraintAssembler<TSparse,TDense>/FinalizeConstraintIteration.md)
- [FinalizeConstraintIteration](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/noop_constraint_assembler/NoOpConstraintAssembler/FinalizeConstraintIteration.md)
- [ExecuteConstraintLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_grid/PGrid<TSparse,TDense>/ExecuteConstraintLoop.md)
- [ExecuteMultigridLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteMultigridLoop.md)
- [ExecuteConstraintLoop](../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/ExecuteConstraintLoop.md)