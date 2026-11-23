---
type: C++ Method
title: MappingPressureToStructure
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L359-L526
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/flags/IsNot
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void MappingPressureToStructure(BinBasedFastPointLocator<3>& node_locator)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeDiscontinuousInterpolation](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/ComputeDiscontinuousInterpolation.md)
- [IsNot](../../../../functions/containers/flags/IsNot.md)