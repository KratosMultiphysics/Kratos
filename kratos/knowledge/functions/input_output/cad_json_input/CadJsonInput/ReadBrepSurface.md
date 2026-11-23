---
type: C++ Method
title: ReadBrepSurface
resource: input_output/cad_json_input.h#L201-L266
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
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBoundaryLoops
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadAndAddEmbeddedEdges
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepSurface( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [GetIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/GetIdOrName.md)
- [HasIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/HasIdOrName.md)
- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [ReadBoundaryLoops](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBoundaryLoops.md)
- [ReadAndAddEmbeddedEdges](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadAndAddEmbeddedEdges.md)

# Called by

- [ReadBrepSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces.md)