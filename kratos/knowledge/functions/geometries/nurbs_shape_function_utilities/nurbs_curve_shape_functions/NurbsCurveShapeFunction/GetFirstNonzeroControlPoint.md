---
type: C++ Method
title: GetFirstNonzeroControlPoint
resource: geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions.h#L141-L144
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateQuadraturePointGeometries
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/GlobalSpaceDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetNonzeroControlPointIndices
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`IndexType GetFirstNonzeroControlPoint() const`

# Called by

- [CreateQuadraturePointGeometries](../../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/CreateQuadraturePointGeometries.md)
- [GlobalSpaceDerivatives](../../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/GlobalSpaceDerivatives.md)
- [GetNonzeroControlPointIndices](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/GetNonzeroControlPointIndices.md)
- [ComputeNurbsShapeFunctionValuesAtSpan](../../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan.md)