---
type: C++ Method
title: ConvertIntersectionPointToEdgeRatio
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L827-L836
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertIntersectionPointToEdgeRatio( const Geometry<Node >& rEdge, const array_1d<double,3>& rIntersectionPoint)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)
- [ComputeExtrapolatedGeometryIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections.md)