---
type: C++ Method
title: IterativeSolver
resource: linear_solvers/iterative_solver.h#L94-L102
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IterativeSolver() : mResidualNorm(0) , mIterationsNumber(0) , mBNorm(0) , mpPreconditioner(new TPreconditionerType()) , mTolerance(0) , mMaxIterationsNumber(0)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)