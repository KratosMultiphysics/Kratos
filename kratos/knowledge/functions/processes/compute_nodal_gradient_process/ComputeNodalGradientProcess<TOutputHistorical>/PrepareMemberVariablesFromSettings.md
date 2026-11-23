---
type: C++ Method
title: PrepareMemberVariablesFromSettings
resource: processes/compute_nodal_gradient_process.cpp#L223-L262
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void ComputeNodalGradientProcess<TOutputHistorical>::PrepareMemberVariablesFromSettings(Parameters ThisParameters)`

# Calls

- [RecursivelyValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/RecursivelyValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)