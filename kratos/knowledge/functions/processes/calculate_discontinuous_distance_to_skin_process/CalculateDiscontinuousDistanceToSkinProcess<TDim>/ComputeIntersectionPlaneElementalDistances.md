---
type: C++ Method
title: ComputeIntersectionPlaneElementalDistances
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L460-L506
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/plane_3d/Plane3D/CalculateSignedDistance
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ReplaceZeroDistances
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeIntersectionPlaneElementalDistances( Element& rElement, const PointerVector<GeometricalObject>& rIntersectedObjects, const std::vector<array_1d<double,3>>& rIntersectionPointsCoordinates)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeElementalDistancesFromPlaneApproximation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromPlaneApproximation.md)
- [CalculateSignedDistance](../../../../functions/geometries/plane_3d/Plane3D/CalculateSignedDistance.md)
- [CorrectDistanceOrientation](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CorrectDistanceOrientation.md)
- [ReplaceZeroDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ReplaceZeroDistances.md)

# Called by

- [CalculateElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalDistances.md)
- [CalculateElementalAndEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances.md)
- [ComputeElementalDistancesFromEdgeRatios](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios.md)