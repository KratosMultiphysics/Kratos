---
type: C++ Method
title: GenerateSkinModelPart
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1533-L1672
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateSkinModelPart( ModelPart& mrNewSkinModelPart )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [prod](../../../../functions/includes/amatrix_interface/prod.md)