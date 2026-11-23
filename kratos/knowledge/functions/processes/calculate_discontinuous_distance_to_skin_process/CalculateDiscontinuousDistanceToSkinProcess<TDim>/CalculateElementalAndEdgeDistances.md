---
type: C++ Method
title: CalculateElementalAndEdgeDistances
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L270-L331
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/SetToSplitFlag
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateElementalAndEdgeDistances( Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)
- [ComputeElementalDistancesFromEdgeRatios](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeElementalDistancesFromEdgeRatios.md)
- [ComputeIntersectionPlaneElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances.md)
- [SetToSplitFlag](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/SetToSplitFlag.md)

# Called by

- [CalculateDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateDistances.md)