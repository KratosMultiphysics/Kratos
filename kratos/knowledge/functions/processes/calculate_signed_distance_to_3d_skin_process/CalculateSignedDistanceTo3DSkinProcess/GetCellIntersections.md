---
type: C++ Method
title: GetCellIntersections
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L2332-L2371
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalculateCoordinateNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int GetCellIntersections(OctreeType::cell_type* cell, double* ray, OctreeType::key_type* ray_key, int direction, std::vector<std::pair<double, Element::GeometryType*> >& intersections)`

# Calls

- [CalculateCoordinateNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalculateCoordinateNormalized.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)