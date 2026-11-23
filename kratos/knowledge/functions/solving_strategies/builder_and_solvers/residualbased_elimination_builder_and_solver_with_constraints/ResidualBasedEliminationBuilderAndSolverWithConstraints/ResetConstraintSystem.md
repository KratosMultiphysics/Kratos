---
type: C++ Method
title: ResetConstraintSystem
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints.h#L1960-L1994
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/containers/csr_matrix/CsrMatrix/nnz
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ResetConstraintSystem()`

# Calls

- [nnz](../../../../../functions/containers/csr_matrix/CsrMatrix/nnz.md)

# Called by

- [BuildWithConstraints](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildWithConstraints.md)
- [BuildRHSNoDirichlet](../../../../../functions/solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_with_constraints/ResidualBasedEliminationBuilderAndSolverWithConstraints/BuildRHSNoDirichlet.md)