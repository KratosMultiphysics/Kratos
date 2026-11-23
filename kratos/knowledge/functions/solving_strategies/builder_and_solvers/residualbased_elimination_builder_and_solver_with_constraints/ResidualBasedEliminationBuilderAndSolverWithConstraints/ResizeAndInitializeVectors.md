---
type: C++ Method
title: ResizeAndInitializeVectors
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1214-L1293
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructRelationMatrixStructure
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeConstraintContribution
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ResizeAndInitializeVectors( typename TSchemeType::Pointer pScheme, TSystemMatrixPointerType& pA, TSystemVectorPointerType& pDx, TSystemVectorPointerType& pb, ModelPart& rModelPart ) override`

# Calls

- [GetReshapeMatrixFlag](../../../../../functions/solving_strategies/builder_and_solvers/builder_and_solver/GetReshapeMatrixFlag.md)
- [ConstructRelationMatrixStructure](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ConstructRelationMatrixStructure.md)
- [ComputeConstraintContribution](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeConstraintContribution.md)