---
type: C++ Method
title: SearchInRadius
resource: spatial_containers/geometrical_objects_bins.cpp#L68-L114
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometricalObjectsBins::SearchInRadius( const PointType& rPoint, const double Radius, std::vector<ResultType>& rResults )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance.md)