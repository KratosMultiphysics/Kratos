---
type: C++ Function
title: IsolateRoots
resource: utilities/polynomial_utilities.cpp#L153-L210
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/polynomial_utilities/Differentiate
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/Sturm
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/SignChanges
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void IsolateRoots( std::vector<IntervalType>& rRootIntervals, const PolynomialType& rPolynomial, const IntervalType& rRange )`

# Calls

- [Differentiate](../../../functions/utilities/polynomial_utilities/Differentiate.md)
- [Sturm](../../../functions/utilities/polynomial_utilities/Sturm.md)
- [SignChanges](../../../functions/utilities/polynomial_utilities/SignChanges.md)
- [pop_back](../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE.md)