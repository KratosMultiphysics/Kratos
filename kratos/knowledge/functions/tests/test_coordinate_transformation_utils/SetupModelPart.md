---
type: Python Function
title: SetupModelPart
resource: tests/test_coordinate_transformation_utils.py#L11-L31
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/variable_utils/SetFlag
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def SetupModelPart(cls, mdpa_name, domain_size):`

# Calls

- [SetFlag](../../../functions/utilities/variable_utils/SetFlag.md)
- [GetGeometry](../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)