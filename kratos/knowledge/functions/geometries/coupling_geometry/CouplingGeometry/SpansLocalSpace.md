---
type: C++ Method
title: SpansLocalSpace
resource: geometries/coupling_geometry.h#L318-L374
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/includes/amatrix_interface/norm_2
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/coupling_geometry/CouplingGeometry/MergeSpans
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual void SpansLocalSpace( std::vector<double>& rSpans, IndexType LocalDirectionIndex = 0) const override`

# Calls

- [norm_2](../../../../functions/includes/amatrix_interface/norm_2.md)
- [MergeSpans](../../../../functions/geometries/coupling_geometry/CouplingGeometry/MergeSpans.md)