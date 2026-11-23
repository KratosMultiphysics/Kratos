---
type: C++ Method
title: ReadBrepCurve
resource: input_output/cad_json_input.h#L436-L459
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/HasIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/GetIdOrName
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepCurve( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [HasIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/HasIdOrName.md)
- [GetIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/GetIdOrName.md)

# Called by

- [ReadBrepEdge](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge.md)