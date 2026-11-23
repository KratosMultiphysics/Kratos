---
type: C++ Method
title: Execute
resource: processes/calculate_distance_to_path_process.cpp#L42-L120
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/variable_utils/SetNonHistoricalVariableToZero
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CalculateDistanceToPathProcess<THistorical>::Execute()`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)
- [SetHistoricalVariableToZero](../../../../functions/utilities/variable_utils/SetHistoricalVariableToZero.md)
- [SetNonHistoricalVariableToZero](../../../../functions/utilities/variable_utils/SetNonHistoricalVariableToZero.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [CalculateDistanceByBruteForce](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceByBruteForce.md)