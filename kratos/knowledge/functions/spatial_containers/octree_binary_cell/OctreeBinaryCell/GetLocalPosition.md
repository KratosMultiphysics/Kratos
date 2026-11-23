---
type: C++ Method
title: GetLocalPosition
resource: spatial_containers/octree_binary_cell.h#L327-L357
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetNodeInNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetNodeInNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t GetLocalPosition(key_type* keys)`

# Called by

- [SetNodeInNeighbours](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetNodeInNeighbours.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes.md)
- [SetNodeInNeighbours](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetNodeInNeighbours.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes.md)