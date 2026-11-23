---
type: C++ Method
title: AssembleRowContributionFreeDofs
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h#L1139-L1202
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/Assemble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRelationMatrix
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleWithoutConstraints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline void AssembleRowContributionFreeDofs( TSystemMatrixType& rA, const Matrix& rALocal, const IndexType i, const IndexType i_local, const Element::EquationIdVectorType& EquationId )`

# Called by

- [Assemble](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver/ResidualBasedEliminationBuilderAndSolver/Assemble.md)
- [AssembleRelationMatrix](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleRelationMatrix.md)
- [AssembleWithoutConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/AssembleWithoutConstraints.md)