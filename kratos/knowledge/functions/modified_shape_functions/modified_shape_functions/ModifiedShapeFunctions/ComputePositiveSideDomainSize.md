---
type: C++ Method
title: ComputePositiveSideDomainSize
resource: modified_shape_functions/modified_shape_functions.cpp#L72-L77
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeDomainSizeOnOneSide
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`double ModifiedShapeFunctions::ComputePositiveSideDomainSize() const`

# Calls

- [ComputeDomainSizeOnOneSide](../../../../functions/modified_shape_functions/modified_shape_functions/ModifiedShapeFunctions/ComputeDomainSizeOnOneSide.md)
- [GetPositiveSubdivisions](../../../../functions/utilities/divide_geometry/DivideGeometry<TPointType>/GetPositiveSubdivisions.md)

# Called by

- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/modified_shape_functions/test_triangle_2d_3_modified_shape_functions/KRATOS_TEST_CASE_IN_SUITE.md)