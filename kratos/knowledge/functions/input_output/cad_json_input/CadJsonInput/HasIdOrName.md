---
type: C++ Method
title: HasIdOrName
resource: input_output/cad_json_input.h#L924-L928
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge
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
---

# Signature

`static bool HasIdOrName( const Parameters rParameters)`

# Called by

- [ReadBrepSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface.md)
- [ReadBrepEdge](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdge.md)
- [ReadBrepCurve](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurve.md)
- [ReadBrepEdgeBrepCurveOnSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface.md)
- [ReadCouplingGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry.md)