---
type: C++ Method
title: ReadBrepEdge
resource: input_output/cad_json_input.h#L411-L434
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/HasIdOrName
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurve
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadBrepEdge( const Parameters rParameters, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [HasIdOrName](../../../../functions/input_output/cad_json_input/CadJsonInput/HasIdOrName.md)
- [ReadBrepCurve](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurve.md)
- [ReadBrepEdgeBrepCurveOnSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface.md)
- [ReadCouplingGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry.md)

# Called by

- [ReadBrepCurveOnSurfaces](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurveOnSurfaces.md)