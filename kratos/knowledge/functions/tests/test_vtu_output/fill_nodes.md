---
type: Python Function
title: fill_nodes
resource: tests/test_vtu_output.py#L215-L227
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_vtu_output/TestVtuOutput/setUpClass
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def fill_nodes(model_part: Kratos.ModelPart):`

# Calls

- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [setUpClass](../../../functions/tests/test_vtu_output/TestVtuOutput/setUpClass.md)