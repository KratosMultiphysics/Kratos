---
type: C++ Method
title: GenerateTetrahedraElements
resource: processes/structured_mesh_generator_process.cpp#L266-L277
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/sources/model_part/ModelPart/CreateNewProperties
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CreateCellTetrahedra
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::GenerateTetrahedraElements()`

# Calls

- [CreateNewProperties](../../../../functions/sources/model_part/ModelPart/CreateNewProperties.md)
- [CreateCellTetrahedra](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CreateCellTetrahedra.md)

# Called by

- [Generate3DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh.md)