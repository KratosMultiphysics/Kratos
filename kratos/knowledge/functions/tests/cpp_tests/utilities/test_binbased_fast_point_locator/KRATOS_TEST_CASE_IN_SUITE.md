---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_binbased_fast_point_locator.cpp#L38-L75
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(BinBasedFastPointLocator1, KratosCoreFastSuite)`

# Calls

- [Create2DGeometry](../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create2DGeometry.md)
- [FindPointOnMeshSimplified](../../../../../functions/utilities/binbased_fast_point_locator/FindPointOnMeshSimplified.md)
- [norm_2](../../../../../functions/includes/amatrix_interface/norm_2.md)
- [Create3DGeometry](../../../../../functions/tests/test_utilities/cpp_tests_utilities/Create3DGeometry.md)