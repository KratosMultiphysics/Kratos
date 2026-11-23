---
type: C++ Method
title: ResizeAndInitializeVectors
resource: solving_strategies/builder_and_solvers/p_multigrid/p_multigrid_builder_and_solver.cpp#L690-L745
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void PMultigridBuilderAndSolver<TSparse,TDense>::ResizeAndInitializeVectors(typename Interface::TSchemeType::Pointer pScheme, typename Interface::TSystemMatrixPointerType& rpLhs, typename Interface::TSystemVectorPointerType& rpSolution, typename Interface::TSystemVectorPointerType& rpRhs, ModelPart& rModelPart)`

# Calls

- [GetReshapeMatrixFlag](../../../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag.md)