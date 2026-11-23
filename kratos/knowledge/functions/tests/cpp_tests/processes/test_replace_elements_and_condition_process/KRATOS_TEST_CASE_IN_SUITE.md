---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/processes/test_replace_elements_and_condition_process.cpp#L101-L133
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTriangle2D3N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTetrahedra3D4N
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/GenerateMixedModelPart
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(ReplaceElementsAndConditionsProcess1, KratosCoreFastSuite)`

# Calls

- [CreateTestModelPartTriangle2D3N](../../../../../functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTriangle2D3N.md)
- [CreateTestModelPartTetrahedra3D4N](../../../../../functions/tests/test_utilities/cpp_tests_utilities/CreateTestModelPartTetrahedra3D4N.md)
- [CreateNewProperties](../../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GenerateMixedModelPart](../../../../../functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/GenerateMixedModelPart.md)
- [GetGeometry](../../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)