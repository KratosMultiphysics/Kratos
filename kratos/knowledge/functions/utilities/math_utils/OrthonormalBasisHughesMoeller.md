---
type: C++ Function
title: OrthonormalBasisHughesMoeller
resource: utilities/math_utils.h#L960-L974
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/UnitCrossProduct
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasis
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void OrthonormalBasisHughesMoeller(const T1& c,T2& a,T3& b )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)
- [UnitCrossProduct](../../../functions/utilities/math_utils/UnitCrossProduct.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [OrthonormalBasis](../../../functions/utilities/math_utils/OrthonormalBasis.md)