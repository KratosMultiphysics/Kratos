---
type: C++ Method
title: ReadGeometryModelPart
resource: input_output/cad_json_input.h#L130-L140
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBreps
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadGeometryModelPart( const Parameters rCadJsonParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [ReadBreps](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBreps.md)

# Called by

- [ReadModelPart](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadModelPart.md)