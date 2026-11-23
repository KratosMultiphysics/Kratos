---
type: Python Method
title: GetNodes
resource: tests/test_vtu_output.py#L451-L459
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def GetNodes(model_part: Kratos.ModelPart, use_model_part_nodes: bool, container: 'typing.Union[Kratos.ConditionsArray, Kratos.ElementsArray]') -> 'list[Kratos.Node]':`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)