---
type: C++ Method
title: Clone
resource: includes/condition.h#L237-L246
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/gid_io/GidIOBase/GetData
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`virtual Pointer Clone (IndexType NewId, NodesArrayType const& ThisNodes) const`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetData](../../../../functions/sources/gid_io/GidIOBase/GetData.md)