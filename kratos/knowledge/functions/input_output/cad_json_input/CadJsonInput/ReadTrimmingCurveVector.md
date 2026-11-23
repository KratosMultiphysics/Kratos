---
type: C++ Method
title: ReadTrimmingCurveVector
resource: input_output/cad_json_input.h#L272-L292
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadTrimmingCurve
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static BrepCurveOnSurfaceLoopType ReadTrimmingCurveVector( const Parameters rParameters, typename NurbsSurfaceType::Pointer pNurbsSurface, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [ReadTrimmingCurve](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadTrimmingCurve.md)