---
type: C++ Method
title: InitializeNodalDistances
resource: processes/calculate_distance_to_skin_process.cpp#L105-L122
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/Initialize
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToSkinProcess<TDim>::InitializeNodalDistances()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [Initialize](../../../../functions/processes/calculate_distance_to_skin_process/CalculateDistanceToSkinProcess<TDim>/Initialize.md)