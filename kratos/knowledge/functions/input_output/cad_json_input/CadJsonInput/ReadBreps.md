---
type: C++ Method
title: ReadBreps
resource: input_output/cad_json_input.h#L146-L178
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadGeometryModelPart
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBreps( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [GetIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/GetIdOrName.md)
- [ReadBrepSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurfaces.md)
- [ReadBrepCurveOnSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces.md)
- [ReadPointsOnGeometries](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries.md)

# Called by

- [ReadGeometryModelPart](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadGeometryModelPart.md)