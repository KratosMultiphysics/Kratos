---
type: C++ Function
title: FromRotationMatrix
resource: utilities/quaternion.h#L475-L526
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_quaternion/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/beam_math_utilities/InverseExponentialTransform
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline Quaternion FromRotationMatrix(const TMatrix3x3 & m)`

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_quaternion/KRATOS_TEST_CASE_IN_SUITE.md)
- [InverseExponentialTransform](../../../functions/utilities/beam_math_utilities/InverseExponentialTransform.md)