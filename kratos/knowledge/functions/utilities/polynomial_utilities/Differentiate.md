---
type: C++ Function
title: Differentiate
resource: utilities/polynomial_utilities.cpp#L88-L100
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/polynomial_utilities/Degree
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/IsolateRoots
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`PolynomialType Differentiate(const PolynomialType& rPolynomial)`

# Calls

- [Degree](../../../functions/utilities/polynomial_utilities/Degree.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [IsolateRoots](../../../functions/utilities/polynomial_utilities/IsolateRoots.md)