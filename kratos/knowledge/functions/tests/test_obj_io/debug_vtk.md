---
type: Python Function
title: debug_vtk
resource: tests/test_obj_io.py#L16-L48
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/python_scripts/project/Project/GetModel
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def debug_vtk(model_part):`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetString](../../../functions/sources/kratos_parameters/Parameters/SetString.md)
- [GetModel](../../../functions/python_scripts/project/Project/GetModel.md)