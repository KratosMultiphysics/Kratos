---
type: C++ Function
title: CollectDoFs
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L321-L345
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/MakeLhsTopology
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void CollectDoFs(const TProxy& rEntities, const ProcessInfo& rProcessInfo, typename Interface::TSchemeType& rScheme, LockObject* pLockBegin, [[maybe_unused]] LockObject* pLockEnd, TIndexSet* pRowSetBegin, [[maybe_unused]] TIndexSet* pRowSetEnd)`

# Called by

- [MakeLhsTopology](../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/PMultigridBuilderAndSolver<TSparse,TDense>/Impl/MakeLhsTopology.md)