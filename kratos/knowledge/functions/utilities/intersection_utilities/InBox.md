---
type: C++ Function
title: InBox
resource: utilities/intersection_utilities.h#L800-L810
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/utilities/intersection_utilities/IntersectionUtilities/ComputeLineBoxIntersection
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline int InBox( const array_1d<double,3>& rIntersectionPoint, const array_1d<double,3>& rBoxPoint0, const array_1d<double,3>& rBoxPoint1, const unsigned int Axis)`

# Called by

- [ComputeLineBoxIntersection](../../../functions/utilities/intersection_utilities/IntersectionUtilities/ComputeLineBoxIntersection.md)