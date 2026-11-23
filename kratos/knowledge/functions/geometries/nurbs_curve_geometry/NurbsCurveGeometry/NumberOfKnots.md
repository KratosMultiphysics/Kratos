---
type: C++ Method
title: NumberOfKnots
resource: geometries/nurbs_curve_geometry.h#L206-L209
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/DomainInterval
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/KnotSpanIntervals
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/IsInsideLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/ClosestPointLocalToLocalSpace
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries/TestGeometry/test_nurbs_curve_3d
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`SizeType NumberOfKnots() const`

# Called by

- [DomainInterval](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/DomainInterval.md)
- [KnotSpanIntervals](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/KnotSpanIntervals.md)
- [IsInsideLocalSpace](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/IsInsideLocalSpace.md)
- [ClosestPointLocalToLocalSpace](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/ClosestPointLocalToLocalSpace.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/geometries/test_nurbs_curve/KRATOS_TEST_CASE_IN_SUITE.md)
- [test_nurbs_curve_3d](../../../../functions/tests/test_geometries/TestGeometry/test_nurbs_curve_3d.md)