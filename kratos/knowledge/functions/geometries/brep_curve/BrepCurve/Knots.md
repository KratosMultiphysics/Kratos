---
type: C++ Method
title: Knots
resource: geometries/brep_curve.h#L211-L214
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/tests/test_geometries/TestGeometry/test_nurbs_curve_3d
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`Vector Knots() const`

# Called by

- [KnotRefinement](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/KnotRefinement.md)
- [DegreeElevation](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation.md)
- [GetBoundaryLoopParameters](../../../../functions/input_output/cad_json_output/CadJsonOutput/GetBoundaryLoopParameters.md)
- [test_nurbs_curve_3d](../../../../functions/tests/test_geometries/TestGeometry/test_nurbs_curve_3d.md)