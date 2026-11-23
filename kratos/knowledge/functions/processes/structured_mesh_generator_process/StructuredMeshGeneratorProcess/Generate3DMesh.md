---
type: C++ Method
title: Generate3DMesh
resource: processes/structured_mesh_generator_process.cpp#L190-L199
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetLocalCoordinatesRange
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes3D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::Generate3DMesh()`

# Calls

- [GetLocalCoordinatesRange](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetLocalCoordinatesRange.md)
- [GenerateNodes3D](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes3D.md)
- [GenerateTetrahedraElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements.md)

# Called by

- [Execute](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Execute.md)