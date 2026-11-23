---
type: C++ Method
title: CalcKeyNormalized
resource: spatial_containers/octree_binary.h#L260-L263
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetRayIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/InsertNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/Insert
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetLeavesInBoundingBoxNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/pGetCellNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/SubvidiveUntilSizeNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`key_type CalcKeyNormalized(coordinate_type coordinate) const`

# Called by

- [GetRayIntersections](../../../../functions/processes/apply_ray_casting_process/ApplyRayCastingProcess<TDim>/GetRayIntersections.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersectionsAndNodes.md)
- [GetIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GetIntersections.md)
- [GetIntersectionsAndNodes](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersectionsAndNodes.md)
- [GetIntersections](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GetIntersections.md)
- [InsertNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/InsertNormalized.md)
- [Insert](../../../../functions/spatial_containers/octree_binary/OctreeBinary/Insert.md)
- [GetLeavesInBoundingBoxNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetLeavesInBoundingBoxNormalized.md)
- [pGetCellNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/pGetCellNormalized.md)
- [SubvidiveUntilSizeNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/SubvidiveUntilSizeNormalized.md)
- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)