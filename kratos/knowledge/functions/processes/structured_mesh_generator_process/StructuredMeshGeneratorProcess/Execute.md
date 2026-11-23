---
type: C++ Method
title: Execute
resource: processes/structured_mesh_generator_process.cpp#L79-L125
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  calls:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/sources/kratos_parameters/Parameters/SetString
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::Execute()`

# Calls

- [Generate2DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh.md)
- [Generate3DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh.md)
- [SetString](../../../../functions/sources/kratos_parameters/Parameters/SetString.md)