---
type: C++ Method
title: GenerateTriangularElements
resource: processes/structured_mesh_generator_process.cpp#L248-L264
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetNodeId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::GenerateTriangularElements()`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [GetNodeId](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetNodeId.md)

# Called by

- [Generate2DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh.md)