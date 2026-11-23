---
type: C++ Method
title: CalculateEmbeddedVelocity
resource: processes/calculate_embedded_signed_distance_to_3d_skin_process.h#L186-L218
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateEmbeddedVelocity(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [Execute](../../../../functions/processes/calculate_embedded_signed_distance_to_3d_skin_process/CalculateEmbeddedSignedDistanceTo3DSkinProcess/Execute.md)