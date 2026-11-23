---
type: Python Method
title: __init__
resource: python_scripts/tikz_output_process.py#L24-L81
visibility: private
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
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/communicator/Communicator/MyPID
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [MyPID](../../../../functions/sources/communicator/Communicator/MyPID.md)
- [DeleteDirectoryIfExisting](../../../../functions/python_scripts/kratos_utilities/DeleteDirectoryIfExisting.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)