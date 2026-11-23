---
type: Python Function
title: GenerateVariableListFromInput
resource: python_scripts/kratos_utilities.py#L113-L116
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/unv_output_process/UnvOutputProcess/ExecuteInitialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GenerateVariableListFromInput(param):`

# Calls

- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)

# Called by

- [ExecuteInitialize](../../../functions/python_scripts/gid_output_process/GiDOutputProcess/ExecuteInitialize.md)
- [ExecuteInitialize](../../../functions/python_scripts/unv_output_process/UnvOutputProcess/ExecuteInitialize.md)