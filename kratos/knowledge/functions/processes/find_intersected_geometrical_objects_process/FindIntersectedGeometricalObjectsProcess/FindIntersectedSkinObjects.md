---
type: C++ Method
title: FindIntersectedSkinObjects
resource: processes/find_intersected_geometrical_objects_process.cpp#L91-L122
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void FindIntersectedGeometricalObjectsProcess::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)`

# Calls

- [GetIntersectedLeaves](../../../../functions/spatial_containers/octree_binary/OctreeBinary/GetIntersectedLeaves.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)