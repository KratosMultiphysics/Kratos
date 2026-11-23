---
type: C++ Function
title: SetUpDistanceField
resource: tests/cpp_tests/processes/test_parallel_distance_calculation_process.cpp#L27-L68
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/processes/test_parallel_distance_calculation_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void SetUpDistanceField( ModelPart& rModelPart, std::function<double(Node& rNode)>& rDistanceFunction, std::function<double&(Node& rNode)>& rDistanceGetter)`

# Calls

- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../../functions/tests/cpp_tests/processes/test_parallel_distance_calculation_process/KRATOS_TEST_CASE_IN_SUITE.md)