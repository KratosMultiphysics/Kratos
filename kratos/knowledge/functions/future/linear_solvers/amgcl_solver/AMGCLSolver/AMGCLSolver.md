---
type: C++ Method
title: AMGCLSolver
resource: future/linear_solvers/amgcl_solver.h#L150-L235
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/future/linear_solvers/amgcl_solver/AMGCLSolver/CheckIfSelectedOptionIsAvailable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`AMGCLSolver(Parameters ThisParameters = Parameters(R"({})"))`

# Calls

- [ValidateAndAssignDefaults](../../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [CheckIfSelectedOptionIsAvailable](../../../../../functions/future/linear_solvers/amgcl_solver/AMGCLSolver/CheckIfSelectedOptionIsAvailable.md)
- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)