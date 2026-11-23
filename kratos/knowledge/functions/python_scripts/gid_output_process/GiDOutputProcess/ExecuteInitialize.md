---
type: Python Method
title: ExecuteInitialize
resource: python_scripts/gid_output_process.py#L156-L233
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GenerateVariableListFromInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/GenerateFlagsListFromInput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
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

`def ExecuteInitialize(self):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetModel](../../../../functions/python_scripts/project/Project/GetModel.md)
- [_initialize_cut_output](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_initialize_cut_output.md)
- [_InitializeGiDIO](../../../../functions/python_scripts/gid_output_process/GiDOutputProcess/_InitializeGiDIO.md)
- [GenerateVariableListFromInput](../../../../functions/python_scripts/kratos_utilities/GenerateVariableListFromInput.md)
- [GenerateFlagsListFromInput](../../../../functions/python_scripts/kratos_utilities/GenerateFlagsListFromInput.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)