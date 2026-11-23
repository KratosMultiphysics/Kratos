---
type: C++ Method
title: GetIntersectionsAndNodes
resource: processes/calculate_signed_distance_to_3d_skin_process.h#L2161-L2272
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetLocalPosition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GetIntersectionsAndNodes(double* ray, int direction, std::vector<std::pair<double,Element::GeometryType*> >& intersections, DistanceSpatialContainersConfigure::data_type& rNodesArray)`

# Calls

- [CalcKeyNormalized](../../../../functions/spatial_containers/octree_binary/OctreeBinary/CalcKeyNormalized.md)
- [GetLocalPosition](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetLocalPosition.md)
- [GetKey](../../../../functions/spatial_containers/octree_binary_cell/OctreeBinaryCell/GetKey.md)