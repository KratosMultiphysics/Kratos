---
type: C++ Method
title: ReadBoundaryLoops
resource: input_output/cad_json_input.h#L327-L369
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/kratos_parameters/Parameters/GetString
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static std::tuple<BrepCurveOnSurfaceLoopArrayType, BrepCurveOnSurfaceLoopArrayType> ReadBoundaryLoops( const Parameters rParameters, typename NurbsSurfaceType::Pointer pNurbsSurface, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [GetString](../../../../functions/sources/kratos_parameters/Parameters/GetString.md)

# Called by

- [ReadBrepSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface.md)