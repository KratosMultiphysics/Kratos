---
type: Python Function
title: components_to_functions
resource: python_scripts/set_initial_state_process.py#L16-L43
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/DependsOnSpace
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def components_to_functions(vect):`

# Calls

- [IsNumber](../../../functions/sources/kratos_parameters/Parameters/IsNumber.md)
- [GetDouble](../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [IsString](../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [GetString](../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [DependsOnSpace](../../../functions/utilities/function_parser_utility/DependsOnSpace.md)

# Called by

- [__init__](../../../functions/python_scripts/set_initial_state_process/SetInitialStateProcess/__init__.md)