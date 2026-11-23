---
type: C++ Function
title: Divide
resource: utilities/polynomial_utilities.cpp#L118-L151
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/polynomial_utilities/Degree
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/Multiply
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/containers/pointer_vector_set/PointerVectorSet/pop_back
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/DropLeadingZeros
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/polynomial_utilities/Sturm
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Divide( PolynomialType& rQuotient, PolynomialType& rRemainder, const PolynomialType& rA, const PolynomialType& rB )`

# Calls

- [Degree](../../../functions/utilities/polynomial_utilities/Degree.md)
- [Multiply](../../../functions/utilities/polynomial_utilities/Multiply.md)
- [pop_back](../../../functions/containers/pointer_vector_set/PointerVectorSet/pop_back.md)
- [DropLeadingZeros](../../../functions/utilities/polynomial_utilities/DropLeadingZeros.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_polynomial_utilities/KRATOS_TEST_CASE_IN_SUITE.md)
- [Sturm](../../../functions/utilities/polynomial_utilities/Sturm.md)