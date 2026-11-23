---
type: C++ Method
title: Tessellate
resource: utilities/tessellation_utilities/curve_tessellation.h#L68-L81
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void Tessellate( const GeometryType& rGeometry, const int PolynomialDegree, const NurbsInterval DomainInterval, const std::vector<NurbsInterval>& rKnotSpanIntervals, const double Tolerance)`

# Calls

- [ComputeTessellation](../../../../../functions/utilities/tessellation_utilities/curve_tessellation/CurveTessellation/ComputeTessellation.md)