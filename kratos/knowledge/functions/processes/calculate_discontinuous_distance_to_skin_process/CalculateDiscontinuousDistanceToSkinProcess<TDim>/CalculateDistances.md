---
type: C++ Method
title: CalculateDistances
resource: processes/calculate_discontinuous_distance_to_skin_process.cpp#L181-L200
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDiscontinuousDistanceToSkinProcess<TDim>::CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)`

# Calls

- [CalculateElementalAndEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CalculateElementalAndEdgeDistances.md)
- [CheckAndCorrectEdgeDistances](../../../../functions/processes/calculate_discontinuous_distance_to_skin_process/CalculateDiscontinuousDistanceToSkinProcess<TDim>/CheckAndCorrectEdgeDistances.md)