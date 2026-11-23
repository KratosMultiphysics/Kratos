---
type: C++ Method
title: IsIntersected
resource: spatial_containers/octree_binary.h#L1090-L1096
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`inline bool IsIntersected(typename cell_type::pointer_type rObject, double Tolerance, const double* rLowPoint, const double* rHighPoint)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)