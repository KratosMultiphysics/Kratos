---
type: Python Function
title: AddVariables
resource: python_scripts/auxiliary_solver_utilities.py#L3-L6
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetStringArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def AddVariables(model_part, aux_variable_list):`

# Calls

- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)