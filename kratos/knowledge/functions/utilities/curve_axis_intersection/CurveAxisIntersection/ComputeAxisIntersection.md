---
type: C++ Method
title: ComputeAxisIntersection
resource: utilities/curve_axis_intersection.h#L146-L224
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputePolygon
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/curve_axis_intersection/CurveAxisIntersection/GetSpanIndex
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/curve_axis_intersection/CurveAxisIntersection/BisectionToAxis
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_curve_axis_intersection/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ComputeAxisIntersection( std::vector<double>& rIntersectionParameters, const GeometryType& rGeometry, double Start, double End, const std::vector<double>& rAxis1, const std::vector<double>& rAxis2, const double Tolerance)`

# Calls

- [ComputePolygon](../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputePolygon.md)
- [GetSpanIndex](../../../../functions/utilities/curve_axis_intersection/CurveAxisIntersection/GetSpanIndex.md)
- [BisectionToAxis](../../../../functions/utilities/curve_axis_intersection/CurveAxisIntersection/BisectionToAxis.md)

# Called by

- [SpansLocalSpace](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/SpansLocalSpace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_curve_axis_intersection/KRATOS_TEST_CASE_IN_SUITE.md)