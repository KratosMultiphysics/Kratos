---
type: C++ Method
title: GetIntersections
resource: processes/calculate_signed_distance_to_3d_condition_skin_process.h#L1793-L1841
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetIntersections(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections)`

# Calls

- [CalcKeyNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized.md)