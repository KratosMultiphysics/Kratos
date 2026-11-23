---
type: C++ Method
title: CalculateCoordinates
resource: spatial_containers/octree_binary.h#L219-L223
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateCoordinates(key_type* keys, coordinate_type * ResultCoordinates) const`

# Called by

- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance.md)