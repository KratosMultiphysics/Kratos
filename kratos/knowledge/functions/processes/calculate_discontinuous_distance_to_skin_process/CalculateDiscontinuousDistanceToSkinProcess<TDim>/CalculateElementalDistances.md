---
type: C++ Method
title: CalculateElementalDistances
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L237-L267
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
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/SetToSplitFlag
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateElementalDistances( Element& rElement1, PointerVector<GeometricalObject>& rIntersectedObjects)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ComputeEdgesIntersections](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeEdgesIntersections.md)
- [ComputeIntersectionPlaneElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances.md)
- [SetToSplitFlag](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/SetToSplitFlag.md)