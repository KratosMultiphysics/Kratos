---
type: C++ Function
title: CreateTestModelPartTriangle2D3N
resource: tests/test_utilities/cpp_tests_utilities.cpp#L82-L120
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateTestModelPartTriangle2D3N(ModelPart& rModelPart)`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_replace_elements_and_condition_process/KRATOS_TEST_CASE_IN_SUITE.md)