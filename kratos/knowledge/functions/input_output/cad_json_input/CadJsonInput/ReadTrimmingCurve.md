---
type: C++ Method
title: ReadTrimmingCurve
resource: input_output/cad_json_input.h#L294-L325
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetBool
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/GetInt
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadTrimmingCurveVector
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static typename BrepCurveOnSurfaceType::Pointer ReadTrimmingCurve( const Parameters rParameters, typename NurbsSurfaceType::Pointer pNurbsSurface, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [GetBool](../../../../functions/sources/kratos_parameters/Parameters/GetBool.md)
- [GetInt](../../../../functions/sources/kratos_parameters/Parameters/GetInt.md)

# Called by

- [ReadTrimmingCurveVector](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadTrimmingCurveVector.md)