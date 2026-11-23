---
type: C++ Function
title: MatrixSquareRoot
resource: utilities/math_utils.h#L1842-L1865
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/GaussSeidelEigenSystem
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/BDBtProductOperation
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool MatrixSquareRoot( const TMatrixType1 &rA, TMatrixType2 &rMatrixSquareRoot, const double Tolerance = 1.0e-18, const SizeType MaxIterations = 20 )`

# Calls

- [GaussSeidelEigenSystem](../../../functions/utilities/math_utils/GaussSeidelEigenSystem.md)
- [BDBtProductOperation](../../../functions/utilities/math_utils/BDBtProductOperation.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)