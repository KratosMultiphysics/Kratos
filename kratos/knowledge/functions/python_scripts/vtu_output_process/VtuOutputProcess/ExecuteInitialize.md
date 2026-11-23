---
type: Python Method
title: ExecuteInitialize
resource: python_scripts/vtu_output_process.py#L83-L100
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/vtu_output/VtuOutput/VtuOutput
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/vtu_output/VtuOutput/AddFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def ExecuteInitialize(self) -> None:`

# Calls

- [VtuOutput](../../../../functions/input_output/vtu_output/VtuOutput/VtuOutput.md)
- [AddVariable](../../../../functions/input_output/vtu_output/VtuOutput/AddVariable.md)
- [GetVariable](../../../../functions/python_interface/kratos_globals/KratosGlobalsImpl/GetVariable.md)
- [AddFlag](../../../../functions/input_output/vtu_output/VtuOutput/AddFlag.md)