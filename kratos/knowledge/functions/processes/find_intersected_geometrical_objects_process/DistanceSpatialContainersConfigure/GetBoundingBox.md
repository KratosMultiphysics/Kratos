---
type: C++ Method
title: GetBoundingBox
resource: processes/find_intersected_geometrical_objects_process.h#L192-L215
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/geometries/geometry/PointsNumber
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline void GetBoundingBox( const typename GeometricalObject::Pointer pObject, double* rLowPoint, double* rHighPoint )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [PointsNumber](../../../../functions/geometries/geometry/PointsNumber.md)