---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/processes/test_set_initial_state_process.cpp#L37-L110
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/constitutive_law/pGetInitialState
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ImposingInitialState2D, KratosCoreFastSuite)`

# Calls

- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [pGetInitialState](../../../../../functions/includes/constitutive_law/pGetInitialState.md)