---
type: C++ Method
title: CalculateDistance
resource: processes/calculate_distance_to_path_process.cpp#L126-L196
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/point_object/PointObject<TObject>/pGetObject
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/FastMinimalDistanceOnLine
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToPathProcess<THistorical>::CalculateDistance( ModelPart& rModelPart, std::vector<Geometry<NodeType>::Pointer>& rVectorSegments )`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [pGetObject](../../../../functions/spatial_containers/point_object/PointObject<TObject>/pGetObject.md)
- [FastMinimalDistanceOnLine](../../../../functions/utilities/geometrical_projection_utilities/FastMinimalDistanceOnLine.md)
- [FastMinimalDistanceOnLineWithRadius](../../../../functions/utilities/geometrical_projection_utilities/GeometricalProjectionUtilities/FastMinimalDistanceOnLineWithRadius.md)