---
type: C++ Function
title: Cofactor
resource: utilities/math_utils.h#L158-L183
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/CofactorMatrix
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double Cofactor(const TMatrixType& rMat, IndexType i, IndexType j)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [CofactorMatrix](../../../functions/utilities/math_utils/CofactorMatrix.md)