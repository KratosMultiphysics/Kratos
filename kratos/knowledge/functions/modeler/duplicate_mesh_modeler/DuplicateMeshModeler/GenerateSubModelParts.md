---
type: C++ Method
title: GenerateSubModelParts
resource: modeler/duplicate_mesh_modeler.cpp#L99-L130
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/HasSubModelPart
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void DuplicateMeshModeler::GenerateSubModelParts( ModelPart& rOriginModelPart, ModelPart& rDestinationModelpart)`

# Calls

- [HasSubModelPart](../../../../functions/sources/model_part/ModelPart/HasSubModelPart.md)

# Called by

- [GenerateMesh](../../../../functions/modeler/duplicate_mesh_modeler/DuplicateMeshModeler/GenerateMesh.md)