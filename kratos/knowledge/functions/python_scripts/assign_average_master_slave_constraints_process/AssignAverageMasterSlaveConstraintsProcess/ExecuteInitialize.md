---
type: Python Method
title: ExecuteInitialize
resource: python_scripts/assign_average_master_slave_constraints_process.py#L31-L50
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__AverageMultiPointConstraint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitialize(self): # Store model part names from the settings`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [__AverageMultiPointConstraint](../../../../functions/python_scripts/assign_average_master_slave_constraints_process/AssignAverageMasterSlaveConstraintsProcess/__AverageMultiPointConstraint.md)