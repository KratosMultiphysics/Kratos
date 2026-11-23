---
type: C++ Method
title: CreateDebugOBB3D
resource: processes/find_intersected_geometrical_objects_with_obb_process.cpp#L416-L443
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasDirectIntersection3D
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIntersectedGeometricalObjectsWithOBBProcess::CreateDebugOBB3D( ModelPart& rModelPart, Properties::Pointer pProperties, OrientedBoundingBox<3>& rOrientedBoundingBox )`

# Called by

- [HasIntersection3D](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection3D.md)
- [HasDirectIntersection3D](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasDirectIntersection3D.md)