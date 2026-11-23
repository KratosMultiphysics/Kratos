---
type: C++ Method
title: ComputeElementalDistancesFromPlaneApproximation
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L442-L457
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/plane_3d/Plane3D/CalculateSignedDistance
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::ComputeElementalDistancesFromPlaneApproximation( Element& rElement, Vector& rElementalDistances, const std::vector<array_1d<double,3>>& rPointVector)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [CalculateSignedDistance](../../../../functions/geometries/plane_3d/Plane3D/CalculateSignedDistance.md)

# Called by

- [ComputeIntersectionPlaneElementalDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/ComputeIntersectionPlaneElementalDistances.md)