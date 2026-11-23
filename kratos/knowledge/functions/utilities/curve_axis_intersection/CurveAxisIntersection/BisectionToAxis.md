---
type: C++ Method
title: BisectionToAxis
resource: utilities/curve_axis_intersection.h#L40-L85
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/curve_axis_intersection/CurveAxisIntersection/ComputeAxisIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static double BisectionToAxis( const GeometryType& rCurve, double IntersectionAxis, double Parameter1, double Parameter2, IndexType AxisDirectionIndex, double Tolerance = 1e-6)`

# Called by

- [ComputeAxisIntersection](../../../../functions/utilities/curve_axis_intersection/CurveAxisIntersection/ComputeAxisIntersection.md)