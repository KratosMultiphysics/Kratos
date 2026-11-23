---
type: C++ Method
title: SearchNearestInRadius
resource: spatial_containers/geometrical_objects_bins.cpp#L119-L157
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/IsObjectFound
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`GeometricalObjectsBins::ResultType GeometricalObjectsBins::SearchNearestInRadius( const PointType& rPoint, const double Radius )`

# Calls

- [SetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/SetDistance.md)
- [SearchNearestInCell](../../../../functions/spatial_containers/geometrical_objects_bins/GeometricalObjectsBins/SearchNearestInCell.md)
- [IsObjectFound](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/IsObjectFound.md)
- [GetDistance](../../../../functions/spatial_containers/spatial_search_result/SpatialSearchResult/GetDistance.md)