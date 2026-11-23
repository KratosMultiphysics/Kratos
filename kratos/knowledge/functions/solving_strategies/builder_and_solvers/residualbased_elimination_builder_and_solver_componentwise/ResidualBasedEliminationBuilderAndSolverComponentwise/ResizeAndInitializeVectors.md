---
type: C++ Method
title: ResizeAndInitializeVectors
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h#L380-L464
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ParallelConstructGraph
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ConstructGraph
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme, TSystemMatrixPointerType& pA, TSystemVectorPointerType& pDx, TSystemVectorPointerType& pb, ModelPart& rModelPart ) override`

# Calls

- [GetReshapeMatrixFlag](../../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag.md)
- [ParallelConstructGraph](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ParallelConstructGraph.md)
- [ConstructGraph](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise/ResidualBasedEliminationBuilderAndSolverComponentwise/ConstructGraph.md)