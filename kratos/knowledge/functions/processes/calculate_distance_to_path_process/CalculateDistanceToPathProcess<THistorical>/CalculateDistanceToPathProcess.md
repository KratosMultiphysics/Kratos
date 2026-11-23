---
type: C++ Method
title: CalculateDistanceToPathProcess
resource: processes/calculate_distance_to_path_process.cpp#L28-L36
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/__base_test_calculate_distance_to_path_process
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`CalculateDistanceToPathProcess<THistorical>::CalculateDistanceToPathProcess( Model& rModel, Parameters ThisParameters ) : mrModel(rModel), mThisParameters(ThisParameters)`

# Calls

- [ValidateAndAssignDefaults](../../../../functions/sources/kratos_parameters/Parameters/ValidateAndAssignDefaults.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [__base_test_calculate_distance_to_path_process](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/__base_test_calculate_distance_to_path_process.md)