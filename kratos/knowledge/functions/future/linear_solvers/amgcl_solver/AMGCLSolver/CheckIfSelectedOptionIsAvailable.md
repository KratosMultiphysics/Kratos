---
type: C++ Method
title: CheckIfSelectedOptionIsAvailable
resource: future/linear_solvers/amgcl_solver.h#L635-L649
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/future/linear_solvers/amgcl_solver/AMGCLSolver/AMGCLSolver
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckIfSelectedOptionIsAvailable( const Parameters ThisParameters, const std::string& rOptionName, const std::set<std::string>& rAvailableOptions)`

# Calls

- [GetString](../../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [AMGCLSolver](../../../../../functions/future/linear_solvers/amgcl_solver/AMGCLSolver/AMGCLSolver.md)