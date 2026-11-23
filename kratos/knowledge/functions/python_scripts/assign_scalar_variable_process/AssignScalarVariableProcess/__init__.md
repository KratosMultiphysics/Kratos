---
type: Python Method
title: __init__
resource: python_scripts/assign_scalar_variable_process.py#L31-L94
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/interval_utility/IntervalUtility/IntervalUtility
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/IsNumber
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/function_parser_utility/DependsOnSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __init__(self, Model, settings ):`

# Calls

- [IntervalUtility](../../../../functions/utilities/interval_utility/IntervalUtility/IntervalUtility.md)
- [IsString](../../../../functions/sources/kratos_parameters/Parameters/IsString.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [IsNumber](../../../../functions/sources/kratos_parameters/Parameters/IsNumber.md)
- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [DependsOnSpace](../../../../functions/utilities/function_parser_utility/DependsOnSpace.md)
- [ApplyFunctionToNodesUtility](../../../../functions/utilities/apply_function_to_nodes_utility/ApplyFunctionToNodesUtility.md)