---
type: C++ Method
title: CalculateElementalDistances
resource: processes/calculate_distance_to_skin_process.cpp#L157-L198
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToSkinProcess<TDim>::CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects)`

# Calls

- [CalculateDistanceToNode](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistanceToNode.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)