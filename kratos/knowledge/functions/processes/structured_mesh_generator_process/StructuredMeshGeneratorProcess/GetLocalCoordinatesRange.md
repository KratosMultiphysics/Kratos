---
type: C++ Method
title: GetLocalCoordinatesRange
resource: processes/structured_mesh_generator_process.cpp#L307-L321
visibility: private
generated:
  by: okf-rs/0.7.0
relationships:
  called_by:
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh
    resolved_by: tree-sitter
    confidence: exact
  - target: functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh
    resolved_by: tree-sitter
    confidence: exact
---

# Signature

`void StructuredMeshGeneratorProcess::GetLocalCoordinatesRange(Point& rMinPoint, Point& rMaxPoint)`

# Called by

- [Generate2DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate2DMesh.md)
- [Generate3DMesh](../../../../functions/processes/structured_mesh_generator_process/StructuredMeshGeneratorProcess/Generate3DMesh.md)