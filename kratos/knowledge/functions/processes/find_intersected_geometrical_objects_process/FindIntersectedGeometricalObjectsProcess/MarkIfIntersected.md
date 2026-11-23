---
type: C++ Method
title: MarkIfIntersected
resource: processes/find_intersected_geometrical_objects_process.cpp#L330-L344
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/IdentifyNearEntitiesAndCheckEntityForIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIntersectedGeometricalObjectsProcess::MarkIfIntersected( GeometricalObject& rIntersectedElement, OtreeCellVectorType& rLeaves )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [IdentifyNearEntitiesAndCheckEntityForIntersection](../../../../functions/processes/find_intersected_geometrical_objects_process/FindIntersectedGeometricalObjectsProcess/IdentifyNearEntitiesAndCheckEntityForIntersection.md)