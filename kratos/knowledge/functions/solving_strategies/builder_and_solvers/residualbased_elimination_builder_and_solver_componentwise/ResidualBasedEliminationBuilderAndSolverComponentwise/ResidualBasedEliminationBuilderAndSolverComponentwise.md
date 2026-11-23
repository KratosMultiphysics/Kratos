---
type: C++ Method
title: ResidualBasedEliminationBuilderAndSolverComponentwise
resource: solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h#L136-L151
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`explicit ResidualBasedEliminationBuilderAndSolverComponentwise( typename TLinearSolver::Pointer pNewLinearSystemSolver, Parameters ThisParameters ) : ResidualBasedEliminationBuilderAndSolverType(pNewLinearSystemSolver)`

# Calls

- [ValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)