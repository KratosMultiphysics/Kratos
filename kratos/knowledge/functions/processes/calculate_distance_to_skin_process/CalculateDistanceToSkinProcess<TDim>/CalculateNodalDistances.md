---
type: C++ Method
title: CalculateNodalDistances
resource: processes/calculate_distance_to_skin_process.cpp#L239-L257
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistances
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToSkinProcess<TDim>::CalculateNodalDistances(NodeScalarGetFunctionType& rGetDistanceFunction)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [CalculateDistances](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/CalculateDistances.md)