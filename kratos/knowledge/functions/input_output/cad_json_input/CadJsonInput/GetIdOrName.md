---
type: C++ Method
title: GetIdOrName
resource: input_output/cad_json_input.h#L909-L921
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBreps
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface
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
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::string GetIdOrName( const Parameters rParameters)`

# Calls

- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)
- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ReadBreps](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBreps.md)
- [ReadBrepSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface.md)
- [ReadBrepCurve](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepCurve.md)
- [ReadBrepEdgeBrepCurveOnSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepEdgeBrepCurveOnSurface.md)
- [ReadCouplingGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadCouplingGeometry.md)
- [ReadPointsOnGeometries](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries.md)