---
type: C++ Method
title: ReadPointsOnGeometries
resource: input_output/cad_json_input.h#L540-L592
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointOnGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBreps
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadPointsOnGeometries( const Parameters rParameters, ModelPart& rModelPart, SizeType EchoLevel = 0)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [GetIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/GetIdOrName.md)
- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [ReadPointOnGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointOnGeometry.md)

# Called by

- [ReadBreps](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBreps.md)