---
type: C++ Method
title: ReconstructSlaveSolutionAfterSolve
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1708-L1747
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/CheckMasterSlaveRelation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ReconstructSlaveSolutionAfterSolve( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemMatrixType& rA, TSystemVectorType& rDx, TSystemVectorType& rb )`

# Calls

- [CheckMasterSlaveRelation](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/CheckMasterSlaveRelation.md)

# Called by

- [BuildAndSolveWithConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildAndSolveWithConstraints.md)