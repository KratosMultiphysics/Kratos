---
type: C++ Function
title: Create2DQuadrilateralsGeometry
resource: tests/test_utilities/cpp_tests_utilities.cpp#L125-L163
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/processes/test_calculate_nodal_area_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_compute_nodal_gradient_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Create2DQuadrilateralsGeometry( ModelPart& rModelPart, const std::string& rEntityName, const bool Initialize, const bool Elements )`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_calculate_nodal_area_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_compute_nodal_gradient_process/KRATOS_TEST_CASE_IN_SUITE.md)