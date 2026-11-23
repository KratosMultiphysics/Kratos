---
type: C++ Method
title: GenerateMesh
resource: modeler/structured_mesh_refinement_modeler.h#L95-L200
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/input_output/cad_json_input/CadJsonInput/GetGeometry
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/spatial_containers/bins_dynamic/BinsDynamic/AddPoint
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateConditions
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void GenerateMesh(ModelPart& rThisModelPart, Element const& rReferenceElement)`

# Calls

- [GetGeometry](../../../../functions/input_output/cad_json_input/CadJsonInput/GetGeometry.md)
- [AddPoint](../../../../functions/spatial_containers/bins_dynamic/BinsDynamic/AddPoint.md)
- [GenerateConditions](../../../../functions/modeler/structured_mesh_refinement_modeler/StructuredMeshRefinementModeler/GenerateConditions.md)