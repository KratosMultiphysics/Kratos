---
type: C++ Method
title: CalculateDistances
resource: processes/calculate_distance_to_skin_process.cpp#L125-L154
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateNodalDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateRayDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToSkinProcess<TDim>::CalculateDistances( std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [CalculateNodalDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateNodalDistances.md)
- [CalculateRayDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateRayDistances.md)