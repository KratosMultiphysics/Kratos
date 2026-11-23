---
type: C++ Method
title: IdentifyNearEntitiesAndCheckEntityForIntersection
resource: processes/find_intersected_geometrical_objects_process.cpp#L318-L325
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/MarkIfIntersected
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIntersectedGeometricalObjectsProcess::IdentifyNearEntitiesAndCheckEntityForIntersection( GeometricalObject::Pointer pGeometricalObject, OtreeCellVectorType& rLeaves )`

# Calls

- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)
- [MarkIfIntersected](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/MarkIfIntersected.md)

# Called by

- [Execute](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/Execute.md)