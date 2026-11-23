---
type: C++ Method
title: MergeSpans
resource: geometries/coupling_geometry.h#L584-L605
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/geometries/coupling_geometry/CouplingGeometry/SpansLocalSpace
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void MergeSpans( std::vector<double>& rResultSpans, const std::vector<double>& rSpans1, const std::vector<double>& rSpans2, double Tolerance = 1e-6)`

# Called by

- [SpansLocalSpace](../../../../functions/geometries/coupling_geometry/CouplingGeometry/SpansLocalSpace.md)