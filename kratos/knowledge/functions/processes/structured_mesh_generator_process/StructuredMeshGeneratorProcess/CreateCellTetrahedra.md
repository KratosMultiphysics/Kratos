---
type: C++ Method
title: CreateCellTetrahedra
resource: processes/structured_mesh_generator_process.cpp#L279-L300
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetNodeId
    resolved_by: tree-sitter
    confidence: exact
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::CreateCellTetrahedra(std::size_t I, std::size_t J, std::size_t K, Properties::Pointer pProperties)`

# Calls

- [GetNodeId](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GetNodeId.md)

# Called by

- [GenerateTetrahedraElements](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/GenerateTetrahedraElements.md)