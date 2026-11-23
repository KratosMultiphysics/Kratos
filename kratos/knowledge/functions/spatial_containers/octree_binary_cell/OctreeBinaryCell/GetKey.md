---
type: C++ Method
title: GetKey
resource: spatial_containers/octree_binary_cell.h#L241-L254
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetNodeInNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetNodeInNeighbours
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMeshNew
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetNeighbourKey
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int GetKey(std::size_t position, key_type* keys) const`

# Called by

- [SetNodeInNeighbours](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/SetNodeInNeighbours.md)
- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes.md)
- [SetNodeInNeighbours](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/SetNodeInNeighbours.md)
- [CalculateNotEmptyLeavesDistance](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/CalculateNotEmptyLeavesDistance.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes.md)
- [PrintGiDMeshNew](../../../../functions/spatial_containers/octree_binary/OctreeBinary/PrintGiDMeshNew.md)
- [GetNeighbourKey](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetNeighbourKey.md)