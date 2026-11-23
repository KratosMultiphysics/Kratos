---
type: C++ Method
title: Generate2DMesh
resource: processes/structured_mesh_generator_process.cpp#L179-L188
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetLocalCoordinatesRange
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes2D
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Execute
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::Generate2DMesh()`

# Calls

- [GetLocalCoordinatesRange](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetLocalCoordinatesRange.md)
- [GenerateNodes2D](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateNodes2D.md)
- [GenerateTriangularElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements.md)

# Called by

- [Execute](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Execute.md)