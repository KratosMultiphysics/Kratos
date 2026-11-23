---
type: C++ Function
title: KRATOS_TEST_CASE_IN_SUITE
resource: tests/cpp_tests/utilities/test_svd_utils.cpp#L28-L54
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/svd_utils/JacobiSingularValueDecomposition
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/svd_utils/SVDConditionNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`KRATOS_TEST_CASE_IN_SUITE(SVDUtilsjacobiSVD2x2Test, KratosCoreFastSuite)`

# Calls

- [JacobiSingularValueDecomposition](../../../../../functions/utilities/svd_utils/JacobiSingularValueDecomposition.md)
- [prod](../../../../../functions/includes/amatrix_interface/prod.md)
- [SVDConditionNumber](../../../../../functions/utilities/svd_utils/SVDConditionNumber.md)