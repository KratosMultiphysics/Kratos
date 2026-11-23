---
type: C++ Method
title: SetBoundingBox
resource: spatial_containers/octree_binary.h#L129-L135
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateOctree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateOctree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/SetOctreeBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/SetOctreeBoundingBox
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetBoundingBox(const coordinate_type * Low, const coordinate_type * High)`

# Called by

- [GenerateOctree](../../../../functions/processes/calculate_signed_distance_to_3d_condition_skin_process/CalculateSignedDistanceTo3DConditionSkinProcess/GenerateOctree.md)
- [GenerateOctree](../../../../functions/processes/calculate_signed_distance_to_3d_skin_process/CalculateSignedDistanceTo3DSkinProcess/GenerateOctree.md)
- [SetOctreeBoundingBox](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/SetOctreeBoundingBox.md)
- [SetOctreeBoundingBox](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/SetOctreeBoundingBox.md)