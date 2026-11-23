---
type: C++ Method
title: AssignMinimalNodalDistance
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L1459-L1479
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistanceFluidStructure
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void AssignMinimalNodalDistance()`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [DistanceFluidStructure](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/DistanceFluidStructure.md)