---
type: C++ Method
title: CheckIfPointIsRepeated
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L426-L439
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
---

# Signature

`bool CalculateDiscontinuousDistanceToSkinProcess<TDim>::CheckIfPointIsRepeated( const array_1d<double,3>& rIntersectionPoint, const std::vector<array_1d<double,3>>& rIntersectionPointsVector, const double& rEdgeTolerance)`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)