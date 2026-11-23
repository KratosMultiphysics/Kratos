---
type: C++ Method
title: CheckIfCutEdgesShareNode
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L847-L879
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`bool CalculateDiscontinuousDistanceToSkinProcess<TDim>::CheckIfCutEdgesShareNode( const Element& rElement, const Element::GeometryType::GeometriesArrayType& rEdgesContainer, const array_1d<double,mNumEdges>& rCutEdgesRatioVector) const`

# Called by

- [ComputeExtrapolatedEdgesIntersectionsIfIncised](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised.md)