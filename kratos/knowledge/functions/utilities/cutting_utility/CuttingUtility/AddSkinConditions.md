---
type: C++ Method
title: AddSkinConditions
resource: utilities/cutting_utility.h#L196-L305
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AddSkinConditions(ModelPart& mr_model_part, ModelPart& mr_new_model_part, int plane_number )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)