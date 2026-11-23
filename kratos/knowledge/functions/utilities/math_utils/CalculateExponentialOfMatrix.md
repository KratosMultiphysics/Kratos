---
type: C++ Function
title: CalculateExponentialOfMatrix
resource: utilities/math_utils.h#L1891-L1916
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/prod
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/Factorial
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/norm_frobenius
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void CalculateExponentialOfMatrix( const TMatrixType& rMatrix, TMatrixType& rExponentialMatrix, const double Tolerance = 1000.0*ZeroTolerance, const SizeType MaxTerms = 200 )`

# Calls

- [prod](../../../functions/includes/amatrix_interface/prod.md)
- [Factorial](../../../functions/utilities/math_utils/Factorial.md)
- [norm_frobenius](../../../functions/includes/amatrix_interface/norm_frobenius.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)