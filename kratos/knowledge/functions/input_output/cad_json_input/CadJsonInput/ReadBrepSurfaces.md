---
type: C++ Method
title: ReadBrepSurfaces
resource: input_output/cad_json_input.h#L184-L199
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/IsArray
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBreps
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepSurfaces( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [IsArray](../../../../functions/sources/kratos_parameters/Parameters/IsArray.md)
- [ReadBrepSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface.md)

# Called by

- [ReadBreps](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBreps.md)