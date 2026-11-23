---
type: C++ Method
title: GenerateElements
resource: modeler/structured_mesh_refinement_modeler.h#L305-L339
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateElements(ModelPart& rThisModelPart, Element& rThisElement, NodesVectorType& VolumeNodesVector, NodesVectorType& SurfaceNodesVector, vector<int>& number_of_divisions, SizeType StartElementId)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)