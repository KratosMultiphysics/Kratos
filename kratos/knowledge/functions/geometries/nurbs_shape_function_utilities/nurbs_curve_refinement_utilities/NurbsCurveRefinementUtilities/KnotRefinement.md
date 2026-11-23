---
type: C++ Method
title: KnotRefinement
resource: geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities.cpp#L17-L195
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/CreateFullKnotVector
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/brep_curve/BrepCurve/Knots
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetUpperSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/CreateReducedKnotVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NurbsCurveRefinementUtilities::KnotRefinement( NurbsCurveGeometryType& rGeometry, std::vector<double>& rKnotsToInsert, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsRefined, Vector& rWeightsRefined)`

# Calls

- [CreateFullKnotVector](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/CreateFullKnotVector.md)
- [Knots](../../../../../functions/geometries/brep_curve/BrepCurve/Knots.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetUpperSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetUpperSpan.md)
- [CreateReducedKnotVector](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/CreateReducedKnotVector.md)