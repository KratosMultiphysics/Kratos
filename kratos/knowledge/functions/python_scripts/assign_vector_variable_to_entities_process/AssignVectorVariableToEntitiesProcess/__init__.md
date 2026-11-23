---
type: Python Method
title: __init__
resource: python_scripts/assign_vector_variable_to_entities_process.py#L21-L79
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/AddEmptyValue
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [PrettyPrintJsonString](../../../../functions/sources/kratos_parameters/Parameters/PrettyPrintJsonString.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [AddValue](../../../../functions/sources/kratos_parameters/Parameters/AddValue.md)
- [AddEmptyValue](../../../../functions/sources/kratos_parameters/Parameters/AddEmptyValue.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)