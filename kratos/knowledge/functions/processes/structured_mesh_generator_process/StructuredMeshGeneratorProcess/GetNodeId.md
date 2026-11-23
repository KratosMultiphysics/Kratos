---
type: C++ Method
title: GetNodeId
resource: processes/structured_mesh_generator_process.cpp#L302-L305
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CreateCellTetrahedra
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`std::size_t StructuredMeshGeneratorProcess::GetNodeId(std::size_t I, std::size_t J, std::size_t K)`

# Called by

- [GenerateTriangularElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTriangularElements.md)
- [CreateCellTetrahedra](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/CreateCellTetrahedra.md)