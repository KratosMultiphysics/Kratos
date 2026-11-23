---
type: C++ Method
title: ComputeExtrapolatedEdgesIntersectionsIfIncised
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L713-L734
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfCutEdgesShareNode
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeExtrapolatedEdgesIntersectionsIfIncised( const Element& rElement, const Element::GeometryType::GeometriesArrayType& rEdgesContainer, unsigned int &rNumCutEdges, array_1d<double,mNumEdges>& rCutEdgesRatioVector, array_1d<double,3> &rExtraGeomNormal, array_1d<double,mNumEdges>& rCutExtraEdgesRatioVector)`

# Calls

- [ComputeExtrapolatedGeometryIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedGeometryIntersections.md)
- [CheckIfCutEdgesShareNode](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfCutEdgesShareNode.md)

# Called by

- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)