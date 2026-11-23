---
type: C++ Method
title: ComputeExtrapolatedGeometryIntersections
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L737-L770
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertEdgeRatioToIntersectionPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputePlaneLineIntersection
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeExtrapolatedGeometryIntersections( const Element& rElement, const Element::GeometryType::GeometriesArrayType& rEdgesContainer, unsigned int& rNumCutEdges, array_1d<double,mNumEdges>& rCutEdgesRatioVector, array_1d<double,3>& rExtraGeomNormal, array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector)`

# Calls

- [ConvertEdgeRatioToIntersectionPoint](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertEdgeRatioToIntersectionPoint.md)
- [ComputePlaneLineIntersection](../../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputePlaneLineIntersection.md)
- [ConvertIntersectionPointToEdgeRatio](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio.md)

# Called by

- [ComputeExtrapolatedEdgesIntersectionsIfIncised](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised.md)