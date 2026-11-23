---
type: C++ Method
title: ComputeBoundaryTessellation
resource: modeler/cad_tessellation_modeler.cpp#L132-L158
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetDouble
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/KnotSpanIntervals
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::vector<array_1d<double, 2>> CadTessellationModeler::ComputeBoundaryTessellation(const BrepCurveOnSurfaceType& rBoundarySegment)`

# Calls

- [GetDouble](../../../../functions/sources/kratos_parameters/Parameters/GetDouble.md)
- [ComputeTessellation](../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation.md)
- [KnotSpanIntervals](../../../../functions/geometries/nurbs_curve_geometry/NurbsCurveGeometry/KnotSpanIntervals.md)

# Called by

- [SetupModelPart](../../../../functions/modeler/cad_tessellation_modeler/CadTessellationModeler/SetupModelPart.md)