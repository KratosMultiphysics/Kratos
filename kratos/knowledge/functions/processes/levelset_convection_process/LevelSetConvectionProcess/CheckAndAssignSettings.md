---
type: C++ Method
title: CheckAndAssignSettings
resource: processes/levelset_convection_process.h#L887-L928
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementRegisteredName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CheckAndAssignSettings(const Parameters ThisParameters)`

# Calls

- [GetConvectionElementRegisteredName](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementRegisteredName.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)