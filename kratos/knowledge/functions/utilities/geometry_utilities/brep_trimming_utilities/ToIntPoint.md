---
type: C++ Function
title: ToIntPoint
resource: utilities/geometry_utilities/brep_trimming_utilities.h#L325-L339
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static Clipper2Lib::Point64 ToIntPoint( const double x, const double y, const double factor)`

# Called by

- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)
- [ComputeSpanTriangulation](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation.md)