---
type: C++ Method
title: ReadAndAddEmbeddedEdges
resource: input_output/cad_json_input.h#L371-L387
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/geometries/brep_surface/BrepSurface/AddEmbeddedEdges
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`static void ReadAndAddEmbeddedEdges( typename BrepSurfaceType::Pointer pBrepSurface, const Parameters rParameters, typename NurbsSurfaceType::Pointer pNurbsSurface, ModelPart& rModelPart, ProjectionAlgorithm ProjectionAlgorithmType, SizeType EchoLevel = 0)`

# Calls

- [AddEmbeddedEdges](../../../../functions/geometries/brep_surface/BrepSurface/AddEmbeddedEdges.md)

# Called by

- [ReadBrepSurface](../../../../functions/input_output/cad_json_input/CadJsonInput/ReadBrepSurface.md)