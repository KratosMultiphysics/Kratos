---
type: C++ Method
title: ComputeSpanTriangulation
resource: utilities/geometry_utilities/brep_trimming_utilities.cpp#L155-L298
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetTessellation
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/ToIntPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/geometries/brep_surface/BrepSurface/ComputeSpanTriangulationLocalSpace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void BrepTrimmingUtilities<TShiftedBoundary>::ComputeSpanTriangulation( const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rOuterLoops, const DenseVector<DenseVector<BrepCurveOnSurfacePointerType>>& rInnerLoops, const double u0, const double u1, const double v0, const double v1, bool& rIsTrimmed, std::vector<Matrix>& rTriangles)`

# Calls

- [GetTessellation](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/GetTessellation.md)
- [ToIntPoint](../../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/ToIntPoint.md)
- [Triangulate_OPT](../../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/Triangulate_OPT.md)

# Called by

- [ComputeSpanTriangulationLocalSpace](../../../../../functions/geometries/brep_surface/BrepSurface/ComputeSpanTriangulationLocalSpace.md)