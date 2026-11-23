---
type: C++ Method
title: AssembleRelationMatrix
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L414-L430
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/AssembleRowContributionFreeDofs
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssembleRelationMatrix( TSystemMatrixType& rT, const LocalSystemMatrixType& rTransformationMatrix, const EquationIdVectorType& rSlaveEquationId, const EquationIdVectorType& rMasterEquationId )`

# Calls

- [AssembleRowContributionFreeDofs](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/AssembleRowContributionFreeDofs.md)