---
type: C++ Method
title: IntersectionBox
resource: processes/find_intersected_geometrical_objects_process.h#L240-L247
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool IntersectionBox( const typename GeometricalObject::Pointer pObject, const PointType& rLowPoint, const PointType& rHighPoint )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)