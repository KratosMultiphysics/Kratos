---
type: C++ Function
title: CreateSphereTriangularMesh
resource: tests/test_utilities/cpp_tests_utilities.cpp#L410-L637
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/node/X0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Y0
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/node/Z0
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/processes/test_find_conditions_neighbours_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void CreateSphereTriangularMesh( ModelPart& rModelPart, const std::string& rConditionName, const double Radius, const std::array<double, 3>& rCenter )`

# Calls

- [X0](../../../../functions/includes/node/X0.md)
- [Y0](../../../../functions/includes/node/Y0.md)
- [Z0](../../../../functions/includes/node/Z0.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_find_conditions_neighbours_process/KRATOS_TEST_CASE_IN_SUITE.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/processes/test_find_nodal_neighbours_process/KRATOS_TEST_CASE_IN_SUITE.md)