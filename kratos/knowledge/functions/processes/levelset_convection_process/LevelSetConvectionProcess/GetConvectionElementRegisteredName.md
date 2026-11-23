---
type: C++ Method
title: GetConvectionElementRegisteredName
resource: processes/levelset_convection_process.h#L930-L948
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementsList
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/levelset_convection_process/LevelSetConvectionProcess/CheckAndAssignSettings
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::string GetConvectionElementRegisteredName(Parameters ThisParameters)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetConvectionElementsList](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementsList.md)
- [GetConvectionElementName](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/GetConvectionElementName.md)

# Called by

- [CheckAndAssignSettings](../../../../functions/processes/levelset_convection_process/LevelSetConvectionProcess/CheckAndAssignSettings.md)