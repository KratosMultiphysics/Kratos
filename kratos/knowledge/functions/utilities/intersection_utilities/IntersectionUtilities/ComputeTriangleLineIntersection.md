---
type: C++ Method
title: ComputeTriangleLineIntersection
resource: utilities/intersection_utilities.cpp#L155-L208
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/math_utils/Norm3
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`int IntersectionUtilities::ComputeTriangleLineIntersection( const array_1d<double,3>& rTrianglePoint1, const array_1d<double,3>& rTrianglePoint2, const array_1d<double,3>& rTrianglePoint3, const array_1d<double,3>& rLinePoint1, const array_1d<double,3>& rLinePoint2, array_1d<double,3>& rIntersectionPoint, const double Epsilon)`

# Calls

- [Norm3](../../../../functions/utilities/math_utils/Norm3.md)