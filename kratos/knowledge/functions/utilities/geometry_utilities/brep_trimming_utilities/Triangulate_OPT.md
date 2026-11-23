---
type: C++ Function
title: Triangulate_OPT
resource: utilities/geometry_utilities/brep_trimming_utilities.h#L89-L245
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/IntPointToDoublePoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/InCone
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/Intersects
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/GetAreaOfTriangle
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void Triangulate_OPT(const Clipper2Lib::Path64& polygon, std::vector<Matrix>& triangles, const double factor, const double span_area)`

# Calls

- [IntPointToDoublePoint](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/IntPointToDoublePoint.md)
- [InCone](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/InCone.md)
- [Intersects](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/Intersects.md)
- [GetAreaOfTriangle](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/GetAreaOfTriangle.md)

# Called by

- [CreateBrepSurfaceTrimmingIntegrationPoints](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/CreateBrepSurfaceTrimmingIntegrationPoints.md)
- [ComputeSpanTriangulation](../../../../functions/utilities/geometry_utilities/brep_trimming_utilities/BrepTrimmingUtilities<TShiftedBoundary>/ComputeSpanTriangulation.md)