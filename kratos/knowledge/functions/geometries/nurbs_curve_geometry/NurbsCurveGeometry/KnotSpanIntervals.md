---
type: C++ Method
title: KnotSpanIntervals
resource: geometries/nurbs_curve_geometry.h#L254-L271
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/cpp_tests/utilities/test_curve_tessellation/KRATOS_TEST_CASE_IN_SUITE
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<NurbsInterval> KnotSpanIntervals() const`

# Calls

- [NumberOfKnots](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/NumberOfKnots.md)

# Called by

- [ComputeBoundaryTessellation](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/ComputeBoundaryTessellation.md)
- [KRATOS_TEST_CASE_IN_SUITE](../../../../functions/tests/cpp_tests/utilities/test_curve_tessellation/KRATOS_TEST_CASE_IN_SUITE.md)