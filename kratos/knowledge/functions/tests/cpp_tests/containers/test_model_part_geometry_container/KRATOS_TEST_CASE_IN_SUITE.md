---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/containers/test_model_part_geometry_container.cpp#L61-L141
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/cpp_tests/containers/test_model_part_geometry_container/GenerateLineModelPartGeometryContainer
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/containers/test_model_part_geometry_container/SetUpTestModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(TestModelPartGeometryContainer, KratosCoreGeometryContainerFastSuite)`

# Calls

- [GenerateLineModelPartGeometryContainer](../../../../../functions/tests/cpp_tests/containers/test_model_part_geometry_container/GenerateLineModelPartGeometryContainer.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [SetUpTestModelPart](../../../../../functions/tests/cpp_tests/containers/test_model_part_geometry_container/SetUpTestModelPart.md)