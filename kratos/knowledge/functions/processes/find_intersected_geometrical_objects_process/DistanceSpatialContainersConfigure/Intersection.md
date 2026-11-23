---
type: C++ Method
title: Intersection
resource: processes/find_intersected_geometrical_objects_process.h#L223-L231
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static inline bool Intersection( const typename GeometricalObject::Pointer pObj1, const typename GeometricalObject::Pointer pObj2 )`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)