---
type: C++ Function
title: OrthonormalBasis
resource: utilities/math_utils.h#L943-L950
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/OrthonormalBasisHughesMoeller
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisFrisvad
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisNaive
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/sources/test_oriented_bounding_box/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentsFromNormal
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void OrthonormalBasis(const T1& c,T2& a,T3& b, const IndexType Type = 0 )`

# Calls

- [OrthonormalBasisHughesMoeller](../../../functions/utilities/math_utils/OrthonormalBasisHughesMoeller.md)
- [OrthonormalBasisFrisvad](../../../functions/utilities/math_utils/OrthonormalBasisFrisvad.md)
- [OrthonormalBasisNaive](../../../functions/utilities/math_utils/OrthonormalBasisNaive.md)

# Called by

- [OrientedBoundingBox](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/sources/test_oriented_bounding_box/KRATOS_TEST_CASE_IN_SUITE.md)
- [ComputeTangentsFromNormal](../../../functions/utilities/mortar_utilities/ComputeTangentsFromNormal.md)
- [ComputeTangentNodeWithSlip](../../../functions/utilities/mortar_utilities/ComputeTangentNodeWithSlip.md)