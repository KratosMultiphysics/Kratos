---
type: C++ Method
title: ReadBrepEdgeBrepCurveOnSurface
resource: input_output/cad_json_input.h#L462-L509
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/HasIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepEdgeBrepCurveOnSurface( const Parameters & rParameters, ModelPart & rModelPart, SizeType EchoLevel = 0)`

# Calls

- [GetIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/GetIdOrName.md)
- [HasIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/HasIdOrName.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)

# Called by

- [ReadBrepEdge](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge.md)