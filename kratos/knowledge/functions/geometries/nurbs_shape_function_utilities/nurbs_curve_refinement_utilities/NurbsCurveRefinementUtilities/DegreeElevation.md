---
type: C++ Method
title: DegreeElevation
resource: geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities.cpp#L197-L543
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_curve/BrepCurve/Knots
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/includes/amatrix_interface/Matrix/fill
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void NurbsCurveRefinementUtilities::DegreeElevation( NurbsCurveGeometryType& rGeometry, SizeType rDegreeToElevate, PointerVector<NodeType>& rPointsRefined, Vector& rKnotsRefined, Vector& rWeightsRefined)`

# Calls

- [Knots](../../../../../functions/geometries/brep_curve/BrepCurve/Knots.md)
- [fill](../../../../../functions/includes/amatrix_interface/Matrix/fill.md)
- [GetBinomCoefficient](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_utilities/GetBinomCoefficient.md)