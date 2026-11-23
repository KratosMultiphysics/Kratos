---
type: C++ Function
title: SVDConditionNumber
resource: utilities/svd_utils.h#L258-L288
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline TDataType SVDConditionNumber( const TInputMatrix& rInputMatrix, const std::string TypeSVD = "Jacobi", const TDataType Tolerance = ZeroTolerance, const IndexType MaxIter = 200 )`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_svd_utils/KRATOS_TEST_CASE_IN_SUITE.md)