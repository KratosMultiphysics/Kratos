---
type: C++ Function
title: UnitCrossProduct
resource: utilities/math_utils.h#L926-L933
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/GetPlanes
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisHughesMoeller
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/math_utils/OrthonormalBasisNaive
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void UnitCrossProduct(T1& c, const T2& a, const T3& b )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [OrientedBoundingBox](../../../functions/geometries/oriented_bounding_box/OrientedBoundingBox<3>/OrientedBoundingBox.md)
- [GetPlanes](../../../functions/geometries/tetrahedra_3d_4/Tetrahedra3D4/GetPlanes.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/utilities/test_math_utils/KRATOS_TEST_CASE_IN_SUITE.md)
- [OrthonormalBasisHughesMoeller](../../../functions/utilities/math_utils/OrthonormalBasisHughesMoeller.md)
- [OrthonormalBasisNaive](../../../functions/utilities/math_utils/OrthonormalBasisNaive.md)