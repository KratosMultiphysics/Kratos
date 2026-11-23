---
type: Python Function
title: AddDofs
resource: python_scripts/auxiliary_solver_utilities.py#L8-L31
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

`def AddDofs(model_part, aux_dofs_list, aux_reaction_list):`

# Calls

- [GetStringArray](../../../functions/sources/kratos_parameters/Parameters/GetStringArray.md)
- [GetVariable](../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)