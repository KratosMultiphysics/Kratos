---
type: C++ Method
title: HasDirectIntersection3D
resource: processes/find_intersected_geometrical_objects_with_obb_process.cpp#L347-L380
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/CreateDebugOBB3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/GetOBBHasIntersectionType
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool FindIntersectedGeometricalObjectsWithOBBProcess::HasDirectIntersection3D( GeometryType& rFirstGeometry, GeometryType& rSecondGeometry )`

# Calls

- [CreateDebugOBB3D](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/CreateDebugOBB3D.md)
- [GetOBBHasIntersectionType](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/GetOBBHasIntersectionType.md)

# Called by

- [HasIntersection](../../../../functions/processes/find_intersected_geometrical_objects_with_obb_process/FindIntersectedGeometricalObjectsWithOBBProcess/HasIntersection.md)