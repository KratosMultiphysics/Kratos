---
type: Python Function
title: CreateFastestAvailableDirectLinearSolver
resource: python_scripts/python_linear_solver_factory.py#L23-L37
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spaces/ublas_space/FastestDirectSolverList
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/eigen_solver_factory/ConstructSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def CreateFastestAvailableDirectLinearSolver(): # Using a default linear solver (selecting the fastest one available)`

# Calls

- [CheckIfApplicationsAvailable](../../../functions/python_scripts/kratos_utilities/CheckIfApplicationsAvailable.md)
- [FastestDirectSolverList](../../../functions/spaces/ublas_space/FastestDirectSolverList.md)

# Called by

- [ConstructSolver](../../../functions/python_scripts/eigen_solver_factory/ConstructSolver.md)