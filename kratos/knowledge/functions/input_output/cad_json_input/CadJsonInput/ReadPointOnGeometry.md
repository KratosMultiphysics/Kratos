---
type: C++ Method
title: ReadPointOnGeometry
resource: input_output/cad_json_input.h#L594-L614
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static GeometryPointerType ReadPointOnGeometry( const Parameters rParameters, ModelPart& rModelPart, GeometryPointerType pBackgroundGeometry, SizeType EchoLevel = 0)`

# Called by

- [ReadPointsOnGeometries](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadPointsOnGeometries.md)