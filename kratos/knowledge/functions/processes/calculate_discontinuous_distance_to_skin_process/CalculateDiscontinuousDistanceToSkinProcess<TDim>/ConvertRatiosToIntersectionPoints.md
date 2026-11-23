---
type: C++ Method
title: ConvertRatiosToIntersectionPoints
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L807-L824
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertEdgeRatioToIntersectionPoint
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ConvertRatiosToIntersectionPoints( const Element::GeometryType& rGeometry, const Element::GeometryType::GeometriesArrayType& rEdgesContainer, const array_1d<double,mNumEdges> &rEdgeRatiosVector, std::vector<array_1d <double,3> > &rIntersectionPointsVector)`

# Calls

- [ConvertEdgeRatioToIntersectionPoint](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertEdgeRatioToIntersectionPoint.md)

# Called by

- [ComputeElementalDistancesFromEdgeRatios](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios.md)