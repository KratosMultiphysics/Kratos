---
type: C++ Method
title: ComputeEdgesIntersections
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L334-L423
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfPointIsRepeated
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`unsigned int CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeEdgesIntersections( Element& rElement1, const PointerVector<GeometricalObject>& rIntersectedObjects, const Element::GeometryType::GeometriesArrayType& rEdgesContainer, array_1d<double,mNumEdges> &rCutEdgesRatioVector, array_1d<double,mNumEdges> &rCutExtraEdgesRatioVector, std::vector<array_1d <double,3> > &rIntersectionPointsArray)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CheckIfPointIsRepeated](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckIfPointIsRepeated.md)
- [ConvertIntersectionPointToEdgeRatio](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ConvertIntersectionPointToEdgeRatio.md)
- [ComputeExtrapolatedEdgesIntersectionsIfIncised](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeExtrapolatedEdgesIntersectionsIfIncised.md)

# Called by

- [CalculateElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances.md)
- [CalculateElementalAndEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances.md)