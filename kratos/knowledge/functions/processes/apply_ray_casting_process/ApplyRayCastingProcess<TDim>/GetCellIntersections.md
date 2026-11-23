---
type: C++ Method
title: GetCellIntersections
resource: processes/apply_ray_casting_process.cpp#L388-L433
visibility: private
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

`int ApplyRayCastingProcess<TDim>::GetCellIntersections( OctreeType::cell_type* cell, const double* ray, OctreeType::key_type* ray_key, const unsigned int direction, std::vector<std::pair<double, Element::GeometryType*> > &rIntersections)`

# Calls

- [CalculateCoordinateNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalculateCoordinateNormalized.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)