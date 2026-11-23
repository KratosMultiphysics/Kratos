---
type: C++ Function
title: FastProjectOnLine2D
resource: utilities/geometrical_projection_utilities.h#L282-L308
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/line_2d_2/IsInside
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/line_2d_2/ProjectionPointGlobalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_line_2d_2/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectOnLine2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline double FastProjectOnLine2D( const TGeometryType& rGeometry, const TPointClass1& rPointToProject, TPointClass2& rPointProjected )`

# Calls

- [norm_2](../../../functions/includes/amatrix_interface/norm_2.md)

# Called by

- [IsInside](../../../functions/geometries/line_2d_2/IsInside.md)
- [ProjectionPointGlobalToLocalSpace](../../../functions/geometries/line_2d_2/ProjectionPointGlobalToLocalSpace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../functions/tests/cpp_tests/geometries/test_line_2d_2/KRATOS_TEST_CASE_IN_SUITE.md)
- [TestFastProjectOnLine2D](../../../functions/tests/cpp_tests/utilities/test_geometrical_projection_utilities/TestFastProjectOnLine2D.md)
- [ProjectedIsInside](../../../functions/utilities/geometry_utilities/GeometryUtils/ProjectedIsInside.md)