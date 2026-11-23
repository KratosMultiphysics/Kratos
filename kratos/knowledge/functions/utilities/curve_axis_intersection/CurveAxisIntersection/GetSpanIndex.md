---
type: C++ Method
title: GetSpanIndex
resource: utilities/curve_axis_intersection.h#L87-L125
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

`static void GetSpanIndex( const std::vector<double>& rAxis, IndexType& rSpanIndex, double& rMin, double& rMax, double Parameter, bool Ascending = true)`

# Called by

- [ComputeAxisIntersection](../../../../functions/utilities/curve_axis_intersection/CurveAxisIntersection/ComputeAxisIntersection.md)