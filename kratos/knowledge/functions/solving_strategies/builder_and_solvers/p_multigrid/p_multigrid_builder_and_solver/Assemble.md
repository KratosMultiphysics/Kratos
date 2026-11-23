---
type: C++ Function
title: Assemble
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L411-L506
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Assemble(ModelPart& rModelPart, typename Interface::TSchemeType& rScheme, std::optional<typename Interface::TSystemMatrixType*> pMaybeLhs, std::optional<typename Interface::TSystemVectorType*> pMaybeRhs)`

# Calls

- [visit](../../../../../functions/python/variable_indexing_python/visit.md)