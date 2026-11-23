---
type: Python Method
title: __base_test_calculate_distance_to_path_process
resource: tests/test_calculate_distance_to_path_process.py#L110-L161
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/SetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceToPathProcess
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/DistanceFunction
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/VTKDebug
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/GidDebug
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_brute_force_zero_radius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_brute_force_radius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_zero_radius
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_radius
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`def __base_test_calculate_distance_to_path_process(self, radius = 0.0, brute_force_calculation = False, check_tolerance = 5e-3, debug = False): # Define model`

# Calls

- [SetDouble](../../../../functions/sources/kratos_parameters/Parameters/SetDouble.md)
- [SetBool](../../../../functions/sources/kratos_parameters/Parameters/SetBool.md)
- [CalculateDistanceToPathProcess](../../../../functions/processes/calculate_distance_to_path_process/CalculateDistanceToPathProcess<THistorical>/CalculateDistanceToPathProcess.md)
- [DistanceFunction](../../../../functions/tests/test_calculate_distance_to_path_process/DistanceFunction.md)
- [VTKDebug](../../../../functions/tests/test_calculate_distance_to_path_process/VTKDebug.md)
- [GidDebug](../../../../functions/tests/test_calculate_distance_to_path_process/GidDebug.md)

# Called by

- [test_calculate_distance_to_path_process_brute_force_zero_radius](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_brute_force_zero_radius.md)
- [test_calculate_distance_to_path_process_brute_force_radius](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_brute_force_radius.md)
- [test_calculate_distance_to_path_process_zero_radius](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_zero_radius.md)
- [test_calculate_distance_to_path_process_radius](../../../../functions/tests/test_calculate_distance_to_path_process/TestCalculateDistanceToPathProcess/test_calculate_distance_to_path_process_radius.md)