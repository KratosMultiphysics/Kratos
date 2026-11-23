---
type: C++ Function
title: GetBinomCoefficient
resource: geometries/nurbs_shape_function_utilities/nurbs_utilities.h#L126-L135
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/GlobalSpaceDerivatives
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static constexpr inline SizeType GetBinomCoefficient(const SizeType N, const SizeType K) noexcept`

# Called by

- [GlobalSpaceDerivatives](../../../../functions/geometries/nurbs_curve_on_surface_geometry/NurbsCurveOnSurfaceGeometry/GlobalSpaceDerivatives.md)
- [DegreeElevation](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_refinement_utilities/NurbsCurveRefinementUtilities/DegreeElevation.md)
- [ComputeNurbsShapeFunctionValuesAtSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_curve_shape_functions/NurbsCurveShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan.md)
- [DegreeElevationU](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationU.md)
- [DegreeElevationV](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_refinement_utilities/NurbsSurfaceRefinementUtilities/DegreeElevationV.md)
- [ComputeNurbsShapeFunctionValuesAtSpan](../../../../functions/geometries/nurbs_shape_function_utilities/nurbs_surface_shape_functions/NurbsSurfaceShapeFunction/ComputeNurbsShapeFunctionValuesAtSpan.md)