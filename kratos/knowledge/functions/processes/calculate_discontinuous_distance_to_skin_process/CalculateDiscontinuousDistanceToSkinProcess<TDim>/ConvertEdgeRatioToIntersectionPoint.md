---
type: C++ Method
title: ConvertEdgeRatioToIntersectionPoint
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L839-L844
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertRatiosToIntersectionPoints
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`array_1d<double,3> CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertEdgeRatioToIntersectionPoint( const Geometry<Node >& rEdge, const double& rEdgeRatio)`

# Called by

- [ComputeExtrapolatedGeometryIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections.md)
- [ConvertRatiosToIntersectionPoints](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertRatiosToIntersectionPoints.md)