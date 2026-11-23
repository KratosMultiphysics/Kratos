---
type: C++ Method
title: CalculateCoordinateNormalized
resource: spatial_containers/octree_binary.h#L196-L199
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetCellIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateCoordinateNormalized(const key_type key, coordinate_type& NormalizedCoordinate) const`

# Called by

- [GetCellIntersections](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetCellIntersections.md)
- [GetCellIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetCellIntersections.md)
- [GetCellIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetCellIntersections.md)