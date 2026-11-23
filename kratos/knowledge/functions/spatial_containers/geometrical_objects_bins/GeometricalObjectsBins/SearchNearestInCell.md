---
type: C++ Method
title: SearchNearestInCell
resource: spatial_containers/geometrical_objects_bins.cpp#L294-L310
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometricalObjectsBins::SearchNearestInCell( const CellType& rCell, const PointType& rPoint, ResultType& rResult, const double MaxRadius )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance.md)
- [SetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance.md)

# Called by

- [SearchNearestInRadius](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInRadius.md)