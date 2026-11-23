---
type: C++ Method
title: ApplyDirichletConditions
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L895-L929
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python/variable_indexing_python/visit
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PMultigridBuilderAndSolver<TSparse,TDense>::ApplyDirichletConditions(typename Interface::TSchemeType::Pointer pScheme, ModelPart& rModelPart, typename Interface::TSystemMatrixType& rLhs, [[maybe_unused]] typename Interface::TSystemVectorType& rSolution, typename Interface::TSystemVectorType& rRhs)`

# Calls

- [visit](../../../../../../functions/python/variable_indexing_python/visit.md)