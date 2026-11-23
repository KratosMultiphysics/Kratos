---
type: C++ Method
title: CalculateSignedDistance
resource: geometries/plane_3d.h#L123-L126
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/triangle_3d_3/TriangleTriangleOverlap
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double CalculateSignedDistance(Point const &rPoint)`

# Called by

- [TriangleTriangleOverlap](../../../../functions/geometries/triangle_3d_3/TriangleTriangleOverlap.md)
- [ComputeElementalDistancesFromPlaneApproximation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation.md)
- [ComputeIntersectionPlaneElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances.md)
- [CalculateDistanceToNode](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode.md)