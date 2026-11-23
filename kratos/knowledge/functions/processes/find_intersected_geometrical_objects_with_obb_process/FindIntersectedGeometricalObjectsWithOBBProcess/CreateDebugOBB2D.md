---
type: C++ Method
title: CreateDebugOBB2D
resource: processes/find_intersected_geometrical_objects_with_obb_process.cpp#L385-L411
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasDirectIntersection2D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIntersectedGeometricalObjectsWithOBBProcess::CreateDebugOBB2D( ModelPart& rModelPart, Properties::Pointer pProperties, OrientedBoundingBox<2>& rOrientedBoundingBox )`

# Called by

- [HasIntersection2D](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection2D.md)
- [HasDirectIntersection2D](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasDirectIntersection2D.md)