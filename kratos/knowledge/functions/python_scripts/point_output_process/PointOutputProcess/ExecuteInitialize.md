---
type: Python Method
title: ExecuteInitialize
resource: python_scripts/point_output_process.py#L68-L107
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__CheckVariableIsSolutionStepVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitialize(self): # getting the ModelPart from the Model`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [__CheckVariableIsSolutionStepVariable](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__CheckVariableIsSolutionStepVariable.md)
- [__SearchPoint](../../../../functions/python_scripts/point_output_process/PointOutputProcess/__SearchPoint.md)