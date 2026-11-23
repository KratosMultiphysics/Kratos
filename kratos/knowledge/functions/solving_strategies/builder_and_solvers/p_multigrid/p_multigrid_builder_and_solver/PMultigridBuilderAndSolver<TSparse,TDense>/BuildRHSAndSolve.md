---
type: C++ Method
title: BuildRHSAndSolve
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L835-L864
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/code_location/CodeLocation/CleanFunctionName
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PMultigridBuilderAndSolver<TSparse,TDense>::BuildRHSAndSolve( typename Interface::TSchemeType::Pointer pScheme, ModelPart& rModelPart, typename Interface::TSystemMatrixType& rLhs, typename Interface::TSystemVectorType& rSolution, typename Interface::TSystemVectorType& rRhs)`

# Calls

- [CleanFunctionName](../../../../../../functions/sources/code_location/CodeLocation/CleanFunctionName.md)