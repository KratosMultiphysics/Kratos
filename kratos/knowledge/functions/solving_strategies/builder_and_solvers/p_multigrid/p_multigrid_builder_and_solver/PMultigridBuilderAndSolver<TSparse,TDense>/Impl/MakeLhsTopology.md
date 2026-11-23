---
type: C++ Method
title: MakeLhsTopology
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L348-L398
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/CollectDoFs
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MakeLhsTopology(const typename Interface::TSchemeType::Pointer& rpScheme, typename Interface::TSystemMatrixType& rLhs, ModelPart& rModelPart)`

# Calls

- [CollectDoFs](../../../../../../../functions/solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver/CollectDoFs.md)
- [visit](../../../../../../../functions/python/variable_indexing_python/visit.md)