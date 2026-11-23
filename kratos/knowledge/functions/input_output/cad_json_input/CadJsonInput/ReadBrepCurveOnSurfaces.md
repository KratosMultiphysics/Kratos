---
type: C++ Method
title: ReadBrepCurveOnSurfaces
resource: input_output/cad_json_input.h#L393-L409
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBreps
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepCurveOnSurfaces( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [ReadBrepEdge](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge.md)

# Called by

- [ReadBreps](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBreps.md)