---
type: C++ Method
title: CheckMasterSlaveRelation
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1642-L1698
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeConstraintContribution
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ReconstructSlaveSolutionAfterSolve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CheckMasterSlaveRelation( typename TSchemeType::Pointer pScheme, ModelPart& rModelPart, TSystemVectorType& rDx, TSystemVectorType& rDxSolved )`

# Calls

- [ComputeConstraintContribution](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ComputeConstraintContribution.md)

# Called by

- [ReconstructSlaveSolutionAfterSolve](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/ReconstructSlaveSolutionAfterSolve.md)