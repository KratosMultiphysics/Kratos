---
type: C++ Method
title: CreateVariablesList
resource: solving_strategies/schemes/residual_based_bdf_custom_scheme.h#L563-L600
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/ResidualBasedBDFCustomScheme
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateVariablesList(Parameters ThisParameters)`

# Calls

- [GetInt](../../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetStringArray](../../../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [ResidualBasedBDFCustomScheme](../../../../../functions/solving_strategies/schemes/residual_based_bdf_custom_scheme/ResidualBasedBDFCustomScheme/ResidualBasedBDFCustomScheme.md)