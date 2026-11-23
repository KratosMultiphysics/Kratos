---
type: C++ Method
title: SearchIsInsideInCell
resource: spatial_containers/geometrical_objects_bins.cpp#L315-L331
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
  called_by:
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GeometricalObjectsBins::SearchIsInsideInCell( const CellType& rCell, const PointType& rPoint, ResultType& rResult )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance.md)

# Called by

- [SearchIsInside](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchIsInside.md)